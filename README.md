# SpacePowerWorkshop.jl

A Julia package for simulating a space-based solar power system, combining orbital mechanics from [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl), thermal modeling, and electrical component simulation using the [Dyad](https://help.juliahub.com/dyad/dev) modeling language.

Presented in the [Power Subsystem Design for Small Satellites in Julia](https://juliahub.com/company/resources/power-subsystem-design-for-small-satellites-webinar) webinar by [JuliaHub](https://juliahub.com).

## What is this thing?

This package contains:
- Dyad components for modeling a power system fed by a solar panel on a satellite
- Julia code for simulating the satellite's orbit and the angle of the satellite with respect to the Sun

Additionally, in the `scripts/` directory, you will find:
- `demo.jl`: Shows how to run the models defined here
- `visualization/simple.jl`: Visualizes the satellite's orbit using Makie.jl and GeoMakie.jl
- `visualization/dashboard.jl`: A 3D visualization of the satellite's orbit and power generation dynamics, plus ground track and satellite perspective.  This is the most useful visualization for the workshop.

## Installation

You will need to set up the [Dyad](https://help.juliahub.com/dyad/dev/) product to use this package.  This requires:
- [a JuliaHub account](https://juliahub.com/ui)
- [installing Dyad Studio, a VS Code extension](https://help.juliahub.com/dyad/dev/installation.html)
- having the `DyadRegistry` available and using the JuliaHub package server (`JULIA_PKG_SERVER="https://juliahub.com"`), which Dyad Studio will set up automatically for you

This package is available through the [Dyad Community Registry](https://github.com/DyadLang/DyadCommunityRegistry).
Simply run this code in a Julia REPL:
```julia
using Pkg
Pkg.Registry.add("https://github.com/DyadLang/DyadCommunityRegistry")
```

Then, you can install this package by running:
```julia
using Pkg
Pkg.add("SpacePowerWorkshop")
```

Alternatively, you can clone the repository and install the dependencies manually:
```bash
# Clone the repository
git clone https://github.com/yourusername/SpacePowerWorkshop.jl.git
cd SpacePowerWorkshop.jl

# Activate the package environment
julia --project -e 'using Pkg; Pkg.activate(".")'

# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Components

### Dyad Components

The package includes several pre-built components defined in the Dyad DSL:

- **TempSensor**: Thermal model accounting for solar irradiance and orbital conditions
- **PVCell**: Photovoltaic cell with nonlinear I-V characteristics
- **SolarPanel**: Complete solar panel system with integrated MPPT
- **DCDC_MPPT**: DC-DC converter with Maximum Power Point Tracking
- **SolarPanelSimple**: Simplified solar panel model for quick simulations

## Visualization

Launch the interactive 3D dashboard to visualize satellite orbits and power generation:

```bash
cd scripts/visualization
# make sure this project is instantiated first
julia --project=. dashboard.jl
```

Features:
- Real-time orbital position tracking
- Power generation graphs
- Eclipse detection visualization
- Interactive camera controls

## Project Structure

```
SpacePowerWorkshop.jl/
├── dyad/                    # Dyad DSL component definitions
├── generated/               # Auto-generated Julia code
├── src/                     # Core Julia source code
│   ├── SpacePowerWorkshop.jl
│   ├── main.jl
│   └── analysis/
├── scripts/                 # Demo and visualization scripts
└── test/                    # Test suite
```
