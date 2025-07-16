# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

SpacePowerWorkshop.jl is a Julia package that simulates space-based solar power systems for satellites. It combines orbital mechanics, thermal modeling, and electrical component simulation using the Dyad domain-specific language.

## Development Commands

```bash
# Activate the package environment
julia --project -e 'using Pkg; Pkg.activate(".")'

# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Run main simulation
julia --project src/main.jl

# Run interactive demo
julia --project scripts/demo.jl

# Run 3D visualization
cd scripts/visualization && julia --project=../.. 3dplot.jl
```

## Architecture

### Core Components

1. **Dyad DSL Models** (`/dyad/`): Domain-specific component definitions that get compiled to Julia
   - `TempSensor.dyad`: Thermal model accounting for solar irradiance and orbital conditions
   - `PVCell.dyad`: Photovoltaic cell with nonlinear I-V characteristics
   - `SolarPanel.dyad`: Complete solar panel system with MPPT
   - `DCDC_MPPT.dyad`: DC-DC converter with Maximum Power Point Tracking

2. **Simulation Engine** (`/src/`):
   - `SpacePowerWorkshop.jl`: Main module that exports all components
   - `main.jl`: Primary simulation script that runs full orbital simulations
   - `orbit_analysis.jl`: Orbital propagation and sun vector calculations using SatelliteToolbox

3. **Generated Code** (`/generated/`): Auto-generated Julia code from Dyad models - DO NOT EDIT MANUALLY

### Key Design Patterns

- **Component Hierarchy**: Dyad components compose into larger systems (PVCell → SolarPanel → SpacecraftPower)
- **Time Interpolation**: Orbital data is interpolated for continuous-time differential equation solving
- **Eclipse Detection**: Uses dot product of sun and position vectors to determine illumination
- **Thermal-Electrical Coupling**: Temperature affects PV performance through temperature coefficients

### Critical Dependencies

- `ModelingToolkit`: Symbolic modeling and equation generation
- `DifferentialEquations`: ODE/DAE solving for circuit simulations
- `SatelliteToolbox`: Orbital mechanics and propagation
- `DyadInterface`: Compiles `.dyad` files to Julia code
- `ElectricalComponents`: Provides electrical component models

### Simulation Flow

1. Define orbital parameters (Keplerian elements)
2. Propagate orbit over time period using SGP4 or J2
3. Calculate sun vectors and eclipse conditions
4. Run coupled thermal-electrical simulation
5. Extract power generation, temperature, and battery state
6. Visualize results with Makie

When modifying simulations, ensure orbital data covers the full simulation timespan and handle eclipse transitions smoothly.