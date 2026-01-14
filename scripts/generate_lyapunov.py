"""
Generate shear flow simulations with Lyapunov exponent computation.

Usage:
    # Single worker
    python scripts/generate_lyapunov.py --output-folder my_lyapunov_output

    # Parallel execution across multiple workers
    python scripts/generate_lyapunov.py -ntot 10 -tid 0  # Worker 0 of 10
    python scripts/generate_lyapunov.py -ntot 10 -tid 1  # Worker 1 of 10
"""

from dotenv import load_dotenv
load_dotenv()

import argparse
from itertools import product
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parents[1]))

from src.global_constants import SF_GRID, OUTPUT_PATH
from src.generate_sf_lyapunov import generate_shear_flow_lyapunov
from src.utils import divide_grid


def get_args():
    parser = argparse.ArgumentParser(
        description='Generate shear flow simulations with Lyapunov exponent computation'
    )

    # Multiprocessing arguments
    parser.add_argument('-ntot', type=int, default=1,
                        help="Total number of tasks")
    parser.add_argument('-tid', type=int, default=0,
                        help="Task ID")

    # Script parameters
    parser.add_argument('--output-folder', type=str,
                        default="shearflow2d_lyapunov",
                        help="Output folder name")

    # Lyapunov-specific parameters
    parser.add_argument('--renorm-interval', type=float, default=0.5,
                        help="Time interval between renormalizations (default: 0.5)")
    parser.add_argument('--perturbation-seed', type=int, default=42,
                        help="Random seed for perturbation initial conditions (default: 42)")
    parser.add_argument('--stop-time', type=float, default=20.0,
                        help="Total simulation time (default: 20.0)")

    # Testing parameters (override grid)
    parser.add_argument('--test', action='store_true',
                        help="Run a quick test with small parameters")

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    # Output directory
    dpath = OUTPUT_PATH / args.output_folder
    dpath.mkdir(parents=True, exist_ok=True)

    # Parameter grid
    if args.test:
        # Small test grid: low resolution, single parameter set
        grid = [{
            'resolution': (128, 256),
            'reynolds': 1e4,
            'schmidt': 1.0,
            'width': 1.0,
            'n_shear': 2,
            'n_blobs': 2,
            'init': "default",
        }]
    else:
        # Use standard SF_GRID
        grid = [dict(zip(SF_GRID.keys(), values))
                for values in product(*SF_GRID.values())]

    safety_factor = 4  # Same as regular shear flow

    # Distribute across workers
    grid = divide_grid(grid, args.ntot, args.tid)
    if grid is None:
        exit()

    # Run simulations
    for kwargs in grid:
        generate_shear_flow_lyapunov(
            **kwargs,
            dpath=dpath,
            safety_factor=safety_factor,
            min_dt=1e-8,
            renorm_interval=args.renorm_interval,
            perturbation_seed=args.perturbation_seed,
            stop_sim_time=args.stop_time,
        )
