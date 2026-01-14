# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository generates physics simulation datasets for "The Well" benchmark collection. It produces two types of 2D fluid dynamics simulations using the Dedalus spectral methods framework:

- **Rayleigh-Bénard Convection (RBC)**: Thermal convection between heated plates
- **Shear Flow (SF)**: Turbulent mixing with tracer transport

## Commands

### Installation
```bash
# System dependencies (Ubuntu)
sudo apt install openmpi-bin openmpi-common libopenmpi-dev libfftw3-dev libfftw3-mpi-dev
export CC=mpicc

# Python dependencies
uv sync
```

### Running Simulations
```bash
# Generate shear flow data (single worker)
python scripts/generate.py --pde-name shearflow2d --output-folder my_output

# Generate Rayleigh-Bénard convection data
python scripts/generate.py --pde-name rbc2d --output-folder my_output

# Parallel execution across multiple workers
python scripts/generate.py --pde-name shearflow2d -ntot 10 -tid 0  # Worker 0 of 10
python scripts/generate.py --pde-name shearflow2d -ntot 10 -tid 1  # Worker 1 of 10
```

## Architecture

### Parameter Grid System
The generation process sweeps over parameter grids defined in `src/global_constants.py`:
- `RBC_GRID`: Rayleigh number, Prandtl number, temperature difference, random seeds
- `SF_GRID`: Reynolds number, Schmidt number, shear layer width, number of shear layers/blobs

The `scripts/generate.py` script expands these grids via Cartesian product and distributes them across workers using `src/utils.divide_grid()`.

### Simulation Modules
- `src/generate_rbc.py`: RBC simulation using Chebyshev (vertical) + Fourier (horizontal) bases with no-slip walls
- `src/generate_sf.py`: Shear flow simulation using doubly-periodic Fourier bases
- `src/generate_*_lyapunov.py`: Extended versions that compute Lyapunov exponents via tangent-space integration

### Output
Simulations write HDF5 files to `output/` containing snapshots of:
- RBC: buoyancy, vorticity, pressure, velocity
- SF: tracer, pressure, velocity, vorticity
