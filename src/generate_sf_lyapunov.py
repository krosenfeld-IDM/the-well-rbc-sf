"""
Generate shear flow simulations with Maximum Lyapunov Exponent computation.

Uses tangent-space integration with periodic renormalization to compute the MLE.
"""
import logging
logger = logging.getLogger(__name__)
import numpy as np
import dedalus.public as d3
import h5py

from src.global_constants import filename_sf_lyapunov


def init_standard_shear(x, z, Lx, n_shear=2, n_blobs=2, width=1.0):
    """Initialize shear flow velocity profile."""
    shear = np.zeros((x.shape[0], z.shape[1]), dtype=np.float32)
    velocity = np.zeros((x.shape[0], z.shape[1]), dtype=np.float32)
    z_shear = np.linspace(-1, 1, n_shear, endpoint=False) + 1/n_shear
    for i, z1 in enumerate(z_shear):
        sign = 2 * (i%2) - 1
        zs = n_shear * (z-z1) / 2 / width
        shear += sign * 1/2 * np.tanh(zs/0.1)
        velocity += 0.1 * np.sin(sign*n_blobs*np.pi*x/Lx) * np.exp(-zs**2/0.01)
    shear += 1/2
    return shear, velocity


def generate_shear_flow_lyapunov(
    resolution,
    reynolds, schmidt,
    n_shear, width, n_blobs, init,
    dpath, safety_factor, min_dt,
    renorm_interval=0.5,
    perturbation_seed=42,
    stop_sim_time=20.0,
):
    """
    Generate shear flow simulation with Maximum Lyapunov Exponent computation.

    Uses tangent-space integration with periodic renormalization.

    Parameters
    ----------
    resolution : tuple
        (Nx, Nz) grid resolution
    reynolds : float
        Reynolds number
    schmidt : float
        Schmidt number
    n_shear : int
        Number of shear layers
    width : float
        Width parameter for shear layers
    n_blobs : int
        Number of perturbation blobs in initial condition
    init : str
        Initial condition type
    dpath : Path
        Output directory
    safety_factor : float
        Safety factor for CFL timestep
    min_dt : float
        Minimum allowed timestep
    renorm_interval : float
        Simulation time interval between renormalizations
    perturbation_seed : int
        Random seed for perturbation initial conditions
    stop_sim_time : float
        Total simulation time (default: 20.0)
    """

    # Parameters
    Lx, Lz = 1, 2
    Nx, Nz = resolution
    dealias = 3/2
    timestepper = d3.SBDF3
    initial_dt = 1e-5
    max_dt = 1.0
    dtype = np.float64

    save_name = filename_sf_lyapunov.format(
        resolution[0], resolution[1], reynolds, schmidt, width, n_shear, n_blobs,
    ).replace('.','_')

    # Bases
    coords = d3.CartesianCoordinates('x', 'z')
    dist = d3.Distributor(coords, dtype=dtype)
    xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
    zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(-Lz/2, Lz/2), dealias=dealias)

    # Fields: Base flow
    p = dist.Field(name='p', bases=(xbasis,zbasis))
    s = dist.Field(name='s', bases=(xbasis,zbasis))
    u = dist.VectorField(coords, name='u', bases=(xbasis,zbasis))
    tau_p = dist.Field(name='tau_p')

    # Fields: Perturbations (tangent space)
    dp = dist.Field(name='dp', bases=(xbasis,zbasis))
    du = dist.VectorField(coords, name='du', bases=(xbasis,zbasis))
    dtau_p = dist.Field(name='dtau_p')

    # Substitutions
    nu = 1 / reynolds  # viscosity
    D = nu / schmidt  # tracer diffusivity
    x, z = dist.local_grids(xbasis, zbasis)
    ex, ez = coords.unit_vector_fields(dist)

    # Problem: Combined base flow + perturbation equations
    problem = d3.IVP([u, s, p, tau_p, du, dp, dtau_p], namespace=locals())

    # Base flow equations
    problem.add_equation("dt(u) + grad(p) - nu*lap(u) = - u@grad(u)")
    problem.add_equation("dt(s) - D*lap(s) = - u@grad(s)")
    problem.add_equation("div(u) + tau_p = 0")
    problem.add_equation("integ(p) = 0")  # Pressure gauge

    # Linearized equations for perturbations
    problem.add_equation("dt(du) + grad(dp) - nu*lap(du) = - u@grad(du) - du@grad(u)")
    problem.add_equation("div(du) + dtau_p = 0")
    problem.add_equation("integ(dp) = 0")  # Pressure gauge for perturbation

    # Solver
    solver = problem.build_solver(timestepper)
    solver.stop_sim_time = stop_sim_time

    # Initial conditions: Base flow
    shear, velocity = init_standard_shear(x, z, Lx, n_shear, n_blobs, width)
    u['g'][0] += shear
    u['g'][1] += velocity

    if init == "default":
        pass
    elif init == "sinusoidal":
        for i in range(Nx):
            u['g'][0][i,:] = np.roll(u['g'][0][i,:], int(10*np.sin(4*np.pi*i/Nx)))
    elif init == "velocity_p":
        zs = n_shear * (z-0.0) / 2 / width
        u['g'][1] += 0.1 * np.sin(2*np.pi*(x-0.5)/Lx) * np.exp(-zs**2/0.01)
    elif init == "velocity_m":
        zs = n_shear * (z-0.0) / 2 / width
        u['g'][1] -= 0.1 * np.sin(2*np.pi*(x-0.5)/Lx) * np.exp(-zs**2/0.01)
    elif init == "random_velocity":
        def rand_min_spaced(T, n, min_space):
            assert 0 < T and min_space * n < T
            x = (T - min_space * n) * np.random.rand(n)
            x = np.sort(x)
            dx = np.diff(x)
            dx += min_space
            x = np.cumsum(np.insert(dx, 0, 0)) + x[0]
            return x
        x_pos = rand_min_spaced(1, n_shear*4, min_space=0.1)
        z_pos = rand_min_spaced(2, n_shear*4, min_space=0.2) - 1
        np.random.shuffle(x_pos)
        np.random.shuffle(z_pos)
        for (x1, z1) in zip(x_pos, z_pos):
            xs = (x-x1) / width
            zs = (z-z1) / width
            u['g'][1] += np.exp(-xs**2/0.02 -zs**2/0.02)
    else:
        raise ValueError(f"Unknown initial condition: {init}")

    # Match the tracer to the shear
    s['g'] = u['g'][0]

    # Initial conditions: Perturbation
    np.random.seed(perturbation_seed)
    du.fill_random('g', distribution='normal', scale=1e-6)
    dp['g'] = 0

    # Normalize perturbation to unit L2 norm
    du.change_scales(1)
    pert_norm_sq = d3.integ(du @ du).evaluate()['g'].flat[0]
    pert_norm = np.sqrt(pert_norm_sq)
    if pert_norm > 0:
        du['g'] /= pert_norm

    # Analysis: Snapshots
    logger.info(f"Output path: {dpath/save_name}")
    snapshots = solver.evaluator.add_file_handler(str(dpath/save_name), sim_dt=0.1, max_writes=200)
    snapshots.add_task(s, name='tracer')
    snapshots.add_task(p, name='pressure')
    snapshots.add_task(u, name='shear_velocity')
    snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')

    # CFL
    CFL = d3.CFL(
        solver, initial_dt=initial_dt, cadence=10, safety=0.2/safety_factor, threshold=0.1,
        max_change=1.5, min_change=0.5, max_dt=max_dt, min_dt=min_dt,
    )
    CFL.add_velocity(u)

    # Flow properties
    flow = d3.GlobalFlowProperty(solver, cadence=10)
    flow.add_property((u@ez)**2, name='w2')

    # Lyapunov computation setup
    log_norms = []          # Log of norms at each renormalization
    renorm_times = []       # Times of renormalization
    running_mle = []        # Running MLE estimates
    next_renorm_time = renorm_interval

    lyap_file = dpath / (save_name + "_lyapunov.h5")

    # Main loop
    try:
        logger.info('Starting main loop with Lyapunov computation')
        while solver.proceed:
            timestep = CFL.compute_timestep()
            solver.step(timestep)

            # Lyapunov renormalization
            if solver.sim_time >= next_renorm_time:
                # Compute perturbation L2 norm
                du.change_scales(1)
                pert_norm_sq = d3.integ(du @ du).evaluate()['g'].flat[0]
                pert_norm = np.sqrt(pert_norm_sq)

                if pert_norm > 0:
                    # Record log(norm) before renormalization
                    log_norms.append(np.log(pert_norm))
                    renorm_times.append(solver.sim_time)

                    # Renormalize perturbation to unit norm
                    du['g'] /= pert_norm
                    dp['g'] = 0

                    # Compute running MLE estimate
                    current_mle = sum(log_norms) / solver.sim_time
                    running_mle.append(current_mle)

                    logger.info(
                        'Renorm: Time=%.4f, norm=%.2e, MLE=%.4f',
                        solver.sim_time, pert_norm, current_mle
                    )

                next_renorm_time += renorm_interval

            # Standard logging
            if (solver.iteration-1) % 100 == 0:
                max_w = np.sqrt(flow.max('w2'))
                logger.info('Iteration=%i, Time=%e, dt=%e, max(w)=%f' %(solver.iteration, solver.sim_time, timestep, max_w))

    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        # Save Lyapunov results
        if len(log_norms) > 0:
            final_mle = sum(log_norms) / solver.sim_time

            with h5py.File(lyap_file, 'w') as f:
                f.create_dataset('mle', data=final_mle)
                f.create_dataset('total_sim_time', data=solver.sim_time)
                f.create_dataset('total_iterations', data=solver.iteration)
                f.create_dataset('renorm_times', data=np.array(renorm_times))
                f.create_dataset('log_norms', data=np.array(log_norms))
                f.create_dataset('running_mle', data=np.array(running_mle))

                # Metadata
                f.attrs['reynolds'] = reynolds
                f.attrs['schmidt'] = schmidt
                f.attrs['n_shear'] = n_shear
                f.attrs['width'] = width
                f.attrs['n_blobs'] = n_blobs
                f.attrs['renorm_interval'] = renorm_interval
                f.attrs['perturbation_seed'] = perturbation_seed

            logger.info('Final MLE: %.6f (saved to %s)', final_mle, lyap_file)

        solver.log_stats()
