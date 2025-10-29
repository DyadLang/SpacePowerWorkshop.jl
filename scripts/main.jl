# This version only uses models implemented in Dyad (recommended)

using SpacePowerWorkshop
using ModelingToolkit
using ModelingToolkit: remake
using Plots
using Plots.PlotMeasures

################ 1. Temperature Sensor ################

@mtkbuild temp = TempSensor() 
prob1 = ODEProblem(temp, [], (0, SpacePowerWorkshop.end_time), [])
sol1 = solve(prob1;  saveat=SpacePowerWorkshop.times_adj)

# Temperature
Plots.plot(sol1.t[1:86400], sol1[temp.T][1:86400],
    title="Temperature", xlabel="Day", ylabel="K", label=false,
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    titlefont=font(18), legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))
# Effective irradiance
Plots.plot(sol1.t[1:86400], sol1[temp.G_eff][1:86400],
    title="Effective Irradiance", xlabel="Day", ylabel="Watt / m²", label=false,
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    titlefont=font(18), legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))


################ 2. PV Cell Validation ################

@mtkbuild validation_test = PVCellValidation() 

# T = 75°C
prob2 = ODEProblem(validation_test, [validation_test.cell.ΔT => 48], (0.0, 10); guesses=[validation_test.cell.Im.i => 8.207600054307171])
sol2 = solve(prob2; dtmax=0.001)
plot(sol2[validation_test.vref.V][1:7717], sol2[-validation_test.i.i][1:7717];
    xlabel="V", ylabel="I", label="75°C",
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))

# T = 50°C
sol3 = solve(remake(prob2; p=[validation_test.cell.ΔT => 23]); dtmax=0.001)
plot!(sol3[validation_test.vref.V][1:8595], sol3[-validation_test.i.i][1:8595], label="50°C", linewidth=3)

# T = 25°C
sol4 = solve(remake(prob2; p=[validation_test.cell.ΔT => -2]); dtmax=0.001)
plot!(sol4[validation_test.vref.V][1:9472], sol4[-validation_test.i.i][1:9472], label="25°C", linewidth=3)

################ 3. Solar Panel Simulation ################

@mtkbuild panel = SolarPanel()
prob5 = ODEProblem(panel, [panel.converter.β => -28, panel.converter.stored_energy => 1], (0.0, SpacePowerWorkshop.end_time); guesses=[panel.cell.Im.i => 8.207600054307171, panel.converter.v => 1])
sol5 = solve(prob5; dtmax=0.001)

##### Day 1 plots #####
# Power produced by PV Cell
plot(sol5.t[1:1650], sol5[panel.cell.Im.I * panel.cell.V.v][1:1650],
    title="Power Produced by PV Cell", xlabel="Day", ylabel="Watts", label=false,
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    titlefont=font(18), legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))
# Stored energy in battery
plot(sol5.t[1:1650], sol5[panel.converter.stored_energy][1:1650],
    title="Battery Stored Energy", xlabel="Day", ylabel="Watt-days", label=false,
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    titlefont=font(18), legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))
