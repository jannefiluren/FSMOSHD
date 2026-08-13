# FlexibleSnowModelOSHD.jl

[ci-img]: https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/jannefiluren/FlexibleSnowModelOSHD.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jannefiluren/FlexibleSnowModelOSHD.jl

[runic-img]: https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black
[runic-url]: https://github.com/fredrikekre/Runic.jl

[zenodo-img]: https://img.shields.io/badge/DOI-10.5281/zenodo.21322824-blue.svg
[zenodo-url]: https://doi.org/10.5281/zenodo.21322824

 [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] [![][runic-img]][runic-url] [![DOI][zenodo-img]][zenodo-url]

A Julia implementation of the **Flexible Snow Model (FSM)** for the **Operational Snow Hydrological Service (OSHD)** at SLF. This package provides a comprehensive snow physics model for simulating snow accumulation and melt processes in complex terrain.

## Overview

FlexibleSnowModelOSHD is a multi-layer snow model that simulates:

- **Snow accumulation and ablation** with detailed physics-based processes
- **Multi-layer snow structure** with dynamic layer evolution 
- **Energy balance** including shortwave/longwave radiation, turbulent fluxes, and ground heat transfer
- **Snow hydraulics** with configurable drainage schemes and liquid water retention
- **Forest canopy interactions** including snow interception, unloading, and subcanopy processes
- **Fractional snow cover** using multiple parameterization approaches

The model is designed for operational snow forecasting applications and supports both point-scale and distributed (gridded) simulations across various surface types including open areas, forests, and glaciers.

## Installation

### Prerequisites
- Julia 1.6 or higher
- Required packages are specified in `Project.toml`

### Installation Steps

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl
   cd FlexibleSnowModelOSHD
   ```

2. **Activate the package environment:**
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. **Load the package:**
   ```julia
   using FlexibleSnowModelOSHD
   ```

## Examples

A simulation representing an open site can be run from the terminal by:

```julia
include("script/run_open_station_example.jl")
```

while a corresponding simulation for a forested site can be run by:

```julia
include("script/run_forest_station_example.jl")
```

## Package Structure

```
FlexibleSnowModelOSHD.jl/
├── src/                            # Source code
│   ├── FlexibleSnowModelOSHD.jl    # Main module
│   ├── types.jl                    # Model data structures  
│   ├── parameters.jl               # Physical constants
│   ├── setup.jl                    # Model initialization
│   ├── step.jl                     # Main physics time step
│   ├── snow.jl                     # Snow physics processes
│   ├── radiation.jl                # Radiation calculations
│   ├── thermal.jl                  # Thermal properties
│   └── ...                         # Additional physics modules
├── script/                         # Simulation scripts
└── test/                           # Unit tests and regression tests
```

