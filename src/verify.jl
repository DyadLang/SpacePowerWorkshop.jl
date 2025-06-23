using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using Plots
using Example1
# using BlockComponents
using DifferentialEquations
using SciMLBase
using ElectricalComponents

num_days = 30
include("orbit_analysis.jl")

@mtkmodel TempSensor begin
    @parameters begin
        G = 1361      # solar irradiance at LEO orbit [W/m^2]
        A = 5            # solar panel area [m^2]
        T_ref = 300        # reference temperature [K]

        α = 0.9            # absorptivity of panel
        ϵ = 0.8            # emissivity of panel
        σ = 5.67e-8       # Stefan-Boltzmann constant [W/(m^2*K^4)]

    end
    @variables begin
        G_eff(t), [output=true]    # effective irradiance at solar panel [W/m^2]
        T(t), [output=true]        # panel temperature [K]
        θ(t)    # angle between solar panel normal and Sun direction [rad]
        in_sunlight(t)  # 1 if satellite is in sunlight, 0 otherwise
    end
    @components begin
        ex_theta = TimeVaryingFunction(f=theta_interp)           # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)     # 1 if satellite is in sunlight, 0 otherwise
    end
    @equations begin
        θ ~ ex_theta.output.u
        in_sunlight ~ ex_sunlight.output.u

        T ~ ((α * G_eff)/(ϵ * σ))^(1/4)
        G_eff ~ max(G * cos(θ) * in_sunlight, 0)
    end
end 

########## PVCell Validation ##########


@mtkmodel PVTestLoad begin
    @parameters begin
        over_v(t) = 0
    end
    
    @components begin 
        cell = PVCell_validate()
        vref = VoltageSource()
        res = Resistor(R=2)
        i = CurrentSensor()
        # src = BlockComponents.Ramp(;start_time=0, offset=0, height=35, duration=10)
        ground = Ground()

        # temp = TempSensor()
    end

    @equations begin 
        cell.over_v ~ over_v
        connect(ground.g, cell.n, res.n)
        connect(cell.p, i.n)
        connect(i.p, res.p)
        # vref.V ~ src.y
        vref.V ~ t*35/10

        # cell.G ~ max(1000*sin(2*pi*t/10), 0)
        # cell.T_reading ~ 25*sin(2*pi*t/2) + 323

        # cell.G ~ temp.G_eff
        # cell.G ~ 1361
        # cell.T_reading ~ max(temp.T, 125)
    end
    @continuous_events begin
        ModelingToolkit.SymbolicContinuousCallback(;eqs=[cell.V.v ~ cell.Vocn], affect=ModelingToolkit.ImperativeAffect(modified=(;over_v, i = cell.Im.I), observed=(;v = cell.V.v, Vocn = cell.Vocn)) do m,o,c,i 
            @show o i.t
            return (;over_v = o.v >= o.Vocn ? 1.0 : 0.0, i=0)
        end, reinitializealg=OrdinaryDiffEq.BrownFullBasicInit())
    end
end

@mtkbuild test2 = PVTestLoad() # additional_passes=[ModelingToolkit.IfLifting]
prob2 = ODEProblem(test2, [], (0.0, end_time); guesses=[test2.cell.Im.i => 8.207600054307171])
sol2 = solve(prob2; dtmax=0.001)


# plot(sol2[test2.vref.V], sol2[-test2.i.i], title="I-V")


########## MPPT Simulation ##########

@mtkmodel MPPTCell begin
    @parameters begin
        over_v(t) = 0
    end
    
    @components begin 
        cell = PVCell() # determines the voltage draw
        res = BetaMPPTLoad() # determines the current draw given the voltage draw
        i = CurrentSensor()
        ground = Ground()
        
        temp = TempSensor()
    end

    @equations begin 
        connect(ground.g, cell.n, res.n)
        connect(cell.p, i.n)
        connect(i.p, res.p)
        connect(cell.Vt, res.Vt)

        cell.G ~ temp.G_eff
        cell.T_reading ~ max(temp.T, 125)
    end
end

@mtkbuild cell = MPPTCell()
prob = ODEProblem(cell, [cell.res.β => -28, cell.res.stored_energy => 1], (0.0, num_days); guesses=[cell.cell.Im.i => 8.207600054307171, cell.res.v => 1])
sol = solve(prob; dtmax=0.001)


plot(sol.t[1:1700], sol[cell.res.stored_energy][1:1700], label="stored_energy")
plot(sol.t[1:1700], sol[cell.cell.Im.I * cell.cell.V.v][1:1700], label="solar panel power")
plot!(sol.t[1:1700], sol[cell.cell.Im.I * cell.cell.V.v - cell.res.hotel_load][1:1700], label="Power - Hotel")
plot!(sol.t[1:1700], sol[cell.res.charge_power][1:1700], title="charge power", label = "charge power")


plot(sol, idxs=[cell.res.stored_energy])
plot(sol, idxs=[cell.cell.Im.I * cell.cell.V.v - cell.res.hotel_load], title="Power - Hotel")





plot(sol, idxs=[cell.cell.T_reading])
plot(sol, idxs=[cell.cell.G])
plot(sol, idxs=[cell.temp.G_eff])
plot(sol, idxs=[cell.cell.Im.I])
plot(sol, idxs=[cell.cell.V.v])
plot(sol, idxs=[cell.cell.Im.I * cell.cell.V.v], title="Power")

# Day 1 data
plot(sol.t[1:1016], sol[cell.cell.T_reading][1:1016])
plot(sol.t[1:1016], sol[cell.cell.G][1:1016])
plot(sol.t[1:1016], sol[cell.temp.G_eff][1:1016])
plot(sol.t[1:1016*29], sol[cell.cell.Im.I][1:1016*29])





