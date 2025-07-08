using Example1
using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using Plots
using Plots.PlotMeasures
using DifferentialEquations
using SciMLBase
using ElectricalComponents
using BlockComponents
using DyadData


################ 1. Orbit Analysis ################

num_days = 30
include("orbit_analysis.jl")

@mtkbuild temp = TempSensor2(;theta=DyadDataset(theta, times_adj), sunlight=DyadDataset(sunlight, times_adj)) 
# DyadDataset([theta, times_adj], "times_adj", "theta")
prob1 = ODEProblem(temp, [], (0, end_time), [])
sol1 = solve(prob1;  saveat=times_adj)

# Temperature
plot(sol1.t[1:86400], sol1[temp.T][1:86400],
    title="Temperature", xlabel="Day", ylabel="K", label=false,
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    titlefont=font(18), legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))
# Effective irradiance
plot(sol1.t[1:86400], sol1[temp.G_eff][1:86400],
    title="Effective Irradiance", xlabel="Day", ylabel="Watt / m²", label=false,
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    titlefont=font(18), legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))


################ 2. PV Cell Validation ################

@mtkmodel PVCellValidation begin
    @components begin 
        cell = PVCell_validate()
        vref = VoltageSource()
        i = CurrentSensor()
        ground = Ground()
    end
    @equations begin 
        connect(ground.g, cell.n, vref.n)
        connect(cell.p, i.n)
        connect(i.p, vref.p)
        vref.V ~ t*35/10
    end
end

@mtkbuild validation_test = PVCellValidation() # additional_passes=[ModelingToolkit.IfLifting]

# T = 75°C
prob2 = ODEProblem(validation_test, [validation_test.cell.ΔT => 48], (0.0, 10); guesses=[validation_test.cell.Im.i => 8.207600054307171])
sol2 = solve(prob2; dtmax=0.001)
plot(sol2[validation_test.vref.V][1:7717], sol2[-validation_test.i.i][1:7717];
    xlabel="V", ylabel="I", label="75°C",
    left_margin=5mm, bottom_margin=5mm, linewidth=3,
    legendfont=font(10), guidefont=font(18), tickfont=font(14),
    size=(800,600))

# T = 50°C
prob3 = ODEProblem(validation_test, [validation_test.cell.ΔT => 23], (0.0, 10); guesses=[validation_test.cell.Im.i => 8.207600054307171])
sol3 = solve(prob3; dtmax=0.001)
plot!(sol3[validation_test.vref.V][1:8595], sol3[-validation_test.i.i][1:8595], label="50°C", linewidth=3)

# T = 25°C
prob4 = ODEProblem(validation_test, [validation_test.cell.ΔT => -2], (0.0, 10); guesses=[validation_test.cell.Im.i => 8.207600054307171])
sol4 = solve(prob4; dtmax=0.001)
plot!(sol4[validation_test.vref.V][1:9472], sol4[-validation_test.i.i][1:9472], label="25°C", linewidth=3)


################ 3. Solar Panel Simulation ################

@mtkmodel SolarPanel begin
    @components begin 
        cell = PVCell()         # Exhibits nonlinear I-V relationship based on temperature and irradiance
        converter = DCDC_MPPT() # Adjusts voltage and current draw from PV Cell (integrated MPPT and battery)
        i = CurrentSensor()
        ground = Ground()
        temp = TempSensor()     # Determines temperature and irradiance from orbital analysis
    end
    @equations begin 
        connect(ground.g, cell.n, converter.n)
        connect(cell.p, i.n)
        connect(i.p, converter.p)
        connect(cell.Vt, converter.Vt)

        cell.G ~ temp.G_eff
        # cell.T_reading ~ max(temp.T, 125)
        cell.T_reading ~ temp.T
    end
end

@mtkbuild panel = SolarPanel()
prob5 = ODEProblem(panel, [panel.converter.β => -28, panel.converter.stored_energy => 1], (0.0, num_days); guesses=[panel.cell.Im.i => 8.207600054307171, panel.converter.v => 1])
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
