using SpacePowerWorkshop

# ModelingToolkit lets us manipulate models defined either in Dyad or Julia
using ModelingToolkit
# DyadInterface contains analysis definitions and tools
using DyadInterface

# ## 1. Temperature Sensor
# First, let's solve and plot the temperature sensor model, 
# to assure ourselves that it works.
# We'll use the `@mtkbuild` macro to get a simplified model.

@time @mtkbuild tempsensor = TempSensor()

# Now, we can run a transient analysis to solve the model:
@time res = TransientAnalysis(; 
    model = tempsensor, 
    stop = SpacePowerWorkshop.end_time,    # defined in orbit_analysis.jl
    saveat = SpacePowerWorkshop.times_adj, # defined in orbit_analysis.jl
)

# Now, we can plot the results.  
# This has a minor workaround for CairoMakie because its scimlbase recipe is a bit outdated
# but in general it will work much better once that is done.
using CairoMakie
sol = rebuild_sol(res) # rebuild_sol from DyadInterface
dense_sol = sol(LinRange(0, 1, 10000); idxs = [tempsensor.T, tempsensor.G_eff])
fig, ax1, plt1 = lines(dense_sol.t, dense_sol[tempsensor.T]; axis = (; ylabel = "K", title = "Temperature"))
ax2, plt2 = lines(fig[2, 1], dense_sol.t, dense_sol[tempsensor.G_eff]; axis = (; ylabel = "W/m²", title = "Effective Irradiance"))
linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)
fig

# ## 2. PV Cell Validation

# For this we will build a model of the PV cell system,
# and validate it against the source paper.

@mtkbuild validation_test = PVCellValidation()
# Here, we'll take a slightly different approach to solving the model.
# Let's only construct it once, and use the `remake` function to re-create
# the problem with different parameters.
prob = ODEProblem(validation_test, [validation_test.cell.ΔT => 48], (0.0, 10); guesses=[validation_test.cell.Im.i => 8.207600054307171])
sol = solve(prob; dtmax=0.001)
# Now, let's plot the results:
f, a, p = plot(
    sol; 
    idxs = (validation_test.vref.V, -validation_test.i.i),
    tspan = (0, sol.t[7717]), label = "75°C",
    axis = (; ylabel = "V", xlabel = "I")
)
# We can use remake to re-create the problem with different parameters.
# Note that this is **much** faster than re-creating the component with a
# different parameter value, since we don't have to re-compile the model.
using ModelingToolkit: remake
@time sol2 = solve(remake(prob; p = [validation_test.cell.ΔT => 23]); dtmax=0.001)
plot!(a, sol2; idxs = (validation_test.vref.V, -validation_test.i.i), tspan = (0, sol2.t[8595]), label = "50°C", color = Cycled(2))
f
@time sol3 = solve(remake(prob; p = [validation_test.cell.ΔT => -2]); dtmax=0.001)
plot!(a, sol3; idxs = (validation_test.vref.V, -validation_test.i.i), tspan = (0, sol3.t[9572]), label = "25°C", color = Cycled(3))
axislegend(a; position = :lb)
f

# ## 3. Solar Panel Simulation
# Finally, let's simulate the full solar panel system.
# This is again defined in `SolarPanel.dyad`, and we can 
# construct it with the `@mtkbuild` macro.

@time @mtkbuild panel = SolarPanel()
prob = ODEProblem(panel, [panel.converter.β => -28, panel.converter.stored_energy => 1], (0.0, SpacePowerWorkshop.end_time); guesses=[panel.cell.Im.i => 8.207600054307171, panel.converter.v => 1])
sol = solve(prob; dtmax=0.001)
# Finally we can plot this too:
f, a1, p1 = plot(sol; idxs = panel.cell.Im.I * panel.cell.V.v, tspan = (0, sol.t[1650]))
a2 = Axis(f[1, 1])
p2 = plot!(a2,sol; linestyle = :solid,idxs = panel.converter.stored_energy, tspan = (0, sol.t[1650]), color = Cycled(2))
linkxaxes!(a1, a2)
hidexdecorations!(a2)
a2.yaxisposition = :right
a2.ylabel = Makie.rich("Stored Energy [J]")
a1.ylabel = Makie.rich("Power [W]")
p2.linewidth = 2
p2.linestyle = :solid
axislegend(a1, [p1, p2], ["Power", "Stored Energy"]; position = :rb)
f