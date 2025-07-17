# Visualizations

This folder contains visualizations for the space power workshop.

- `simple.jl` just shows the orbit of the satellite.
- `dashboard.jl` is an animated, interactive dashboard that shows the scenario that was shown in the workshop.  You can't change anything or tune any parameters.
- `dashboard_orbit.jl` TODO allow changing orbit, play button

## How to run

To set up, clone this repo and cd to its top level directory.  Then, run:
```bash
julia --project=./scripts/visualization -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

This will instantiate the environment, and must be run on Julia v1.11 or higher (only because the project uses `[sources]`.).

Then, to run a dashboard, simply run that script:
```bash
julia --startup=no --project=scripts/visualization scripts/visualization/dashboard.jl
```

If you want to interact, run it with `-i`.  Otherwise the Julia session will stop when you close the dashboard.
