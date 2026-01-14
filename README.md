This repository contains the code used to generate the following datasets from "the Well" [1]:
- Rayleigh-Bénard convection
- shear flow

It relies heavily on the dedalus package https://dedalus-project.org/ [2]. 


# Installation

Conda:
```
eval "$(/home/krosenfeld/miniforge3/bin/conda shell.bash hook)"
conda create -n dedalus3
conda activate dedalus3
conda env config vars set OMP_NUM_THREADS=1
conda env config vars set NUMEXPR_MAX_THREADS=1
conda install -c conda-forge dedalus
```

From a venv with python>=3.11.2 run the commands below to install the required packages
```bash
python -m venv ~/venvs/thewell
source ~/venvs/thewell/bin/activate
pip install -r requirements.txt
```

# Generation

The generation is performed by running `generate_rbc.py` or `generate_sf.py`. 
The grid of parameters used (for the PDE and initial conditions) is present in `global_constants.py`.

```python
RBC_GRID = {
    'resolution': [(512, 128)],
    'rayleigh': [1e6, 1e7, 1e8, 1e9, 1e10],
    'prandtl': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    'dT': [0.2, 0.4, 0.6, 0.8, 1.0],
    'seed': [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
    'init': ["default"],
}

SF_GRID = {
    'resolution': [(1024, 2048)],  # downsampled to (256, 512)
    'reynolds': [1e4, 5e4, 1e5, 5e5],
    'schmidt': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    'width': [0.25, 0.5, 1.0, 2.0, 4.0],
    'n_shear': [2, 4],
    'n_blobs': [2, 3, 4, 5],
    'init': ["default"],
}

```

# Lyapunov Exponent Computation

The repository also supports computing the Maximum Lyapunov Exponent (MLE) for shear flow simulations using tangent-space integration.

```bash
# Run with full parameter grid
python scripts/generate_lyapunov.py --output-folder lyapunov_output

# Quick test with reduced resolution (128x256) and short simulation
python scripts/generate_lyapunov.py --test --stop-time 2.0 --output-folder test_lyapunov

# Parallel execution across workers
python scripts/generate_lyapunov.py -ntot 10 -tid 0  # Worker 0 of 10
```

Options:
- `--stop-time`: Total simulation time (default: 20.0)
- `--renorm-interval`: Time between perturbation renormalizations (default: 0.5)
- `--perturbation-seed`: Random seed for perturbation initial conditions (default: 42)
- `--test`: Use small test parameters (128x256 resolution, single parameter set)

Output files:
- `{name}_lyapunov.h5`: Contains MLE value, time series of running estimates, and metadata

# Installation

```
sudo apt update
sudo apt install openmpi-bin openmpi-common libopenmpi-dev
export CC=mpicc
sudo apt install -y libfftw3-dev libfftw3-mpi-dev
uv sync
export UV_ENV_FILE="/home/krosenfeld/projects/the-well-rbc-sf/.env"
```

# References

[1] "The Well: a Large-Scale Collection of Diverse Physics Simulations for Machine Learning"
R.Ohana*, M.Mccabe*, L.Meyer, R.Morel et al. https://arxiv.org/abs/2412.00568

[2] "Dedalus: A flexible framework for numerical simulations with spectral methods" \
K.J.Burns et al. https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.2.023068