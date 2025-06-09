using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using Plots
using Example1
using BlockComponents
using DifferentialEquations
using SciMLBase

end_time = 30
include("orbit_analysis.jl")

# @mtkbuild verify = PVArrayVerification()
# prob = ODEProblem(verify,[verify.pVArray.i => 0],(0,1),[])
# sol = solve(prob)


# @mtkbuild cell = PVCell()
# prob = ODEProblem(cell,[],(0,1),[])
# sol = solve(prob)

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
        ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
        # MPPT 
        # DC-DC converter (takes in voltage)
    end
    @equations begin
        θ ~ ex_theta.output.u
        in_sunlight ~ ex_sunlight.output.u

        T ~ ((α * G_eff)/(ϵ * σ))^(1/4)
        G_eff ~ max(G * cos(θ) * in_sunlight, 0)
    end
end 

@mtkbuild temp = TempSensor() 
prob = ODEProblem(temp, [], (0, end_time), [])
sol = solve(prob;  saveat=times_adj)

plot(sol.t[1:43200], sol[temp.T][1:43200], label="temp.T")
plot(sol.t[1:43200], sol[temp.G_eff][1:43200], label="temp.G_eff")
plot(sol.t[1:43200], sol[temp.θ][1:43200], label="temp.θ")
plot(sol.t[1:43200], sol[temp.in_sunlight][1:43200], label="temp.in_sunlight")
# plot(sol.t[1:43200], sol[temp.sun_facing][1:43200], label="temp.sun_facing")
plot(sol.t[1:43200], sol[cos(temp.θ)][1:43200], label="cos(θ)")

@mtkmodel PVTest begin
    @parameters begin
        over_v(t) = 0
    end
    
    @components begin 
        cell = PVCell()
        vref = VoltageSource()
        i = CurrentSensor()
        # src = BlockComponents.Ramp(;start_time=0, offset=0, height=35, duration=10)
        ground = Ground()

        # temp = TempSensor()
    end

    @equations begin 
        cell.over_v ~ over_v
        connect(ground.g, cell.n, vref.n)
        connect(cell.p, i.n)
        connect(i.p, vref.p)
        # vref.V ~ src.y
        vref.V ~ t*35/10

        cell.G ~ 200
        # cell.G ~ temp.G_eff
        # temp.T ~ cell.T_reading
    end
    @continuous_events begin
        (cell.V.v ~ cell.Vocn) => ModelingToolkit.ImperativeAffect(modified=(;over_v, i = cell.Im.I), observed=(;v = vref.V, Vocn = cell.Vocn)) do m,o,c,i 
            @show o i.t
            return (;over_v = o.v >= o.Vocn ? 1.0 : 0.0, i=0)
        end
    end
end

@mtkbuild test = PVTest() # additional_passes=[ModelingToolkit.IfLifting]
# prob = ODEProblem(test, [], (0.0, end_time))
# prob = ODEProblem(test, [], (0.0, end_time); guesses=[test.cell.Rp_c.i => -0.004364211122503292, test.cell.Im.i => 8.207600054307171])
prob = ODEProblem(test, [], (0.0, end_time); guesses=[test.cell.Im.i => 8.207600054307171])
# sol = solve(prob; dtmax=0.001, initializealg = SciMLBase.OverrideInit() )
sol = solve(prob; dtmax=0.001)

plot(sol[test.vref.V], sol[-test.i.i], title="I-V")

plot(sol, idxs=[test.cell.V.v])
plot(sol, idxs=[test.vref.V])
plot(sol, idxs=[test.src.y])

plot(sol, idxs=[test.cell.Im.I])
plot(sol, idxs=[test.cell.I.i])
plot(sol, idxs=[test.i.i])

plot(sol, idxs=[test.cell.T_reading])
plot(sol, idxs=[test.cell.G])
plot(sol, idxs=[test.temp.T])
plot(sol, idxs=[test.temp.G_eff])
plot(sol, idxs=[test.temp.θ])

plot(sol, idxs=[test.cell.rolloff])
plot(sol, idxs=[test.cell.ipv])



@mtkmodel PVTestLoad begin
    @parameters begin
        over_v(t) = 0
    end
    
    @components begin 
        cell = PVCell()
        # vref = VoltageSource()
        res = Resistor(R=2)
        i = CurrentSensor()
        # src = BlockComponents.Ramp(;start_time=0, offset=0, height=35, duration=10)
        ground = Ground()

        temp = TempSensor()
    end

    @equations begin 
        cell.over_v ~ over_v
        connect(ground.g, cell.n, res.n)
        connect(cell.p, i.n)
        connect(i.p, res.p)
        # vref.V ~ src.y
        # vref.V ~ t*35/10

        # cell.G ~ max(1000*sin(2*pi*t/10), 0)
        # cell.T_reading ~ 25*sin(2*pi*t/2) + 323

        cell.G ~ temp.G_eff
        cell.T_reading ~ max(temp.T, 125)
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
plot(sol2, idxs=[test2.cell.T_reading])
plot(sol2, idxs=[test2.cell.G])
plot(sol2, idxs=[test2.temp.G_eff])
plot(sol2, idxs=[test2.cell.Im.I])

# Day 1 data
plot(sol2.t[1:1016], sol2[test2.cell.T_reading][1:1016])
plot(sol2.t[1:1016], sol2[test2.cell.G][1:1016])
plot(sol2.t[1:1016], sol2[test2.temp.G_eff][1:1016])
plot(sol2.t[1:1016*29], sol2[test2.cell.Im.I][1:1016*29])


plot(sol2, idxs=[test2.cell.V.v])
plot(sol2, idxs=[test2.vref.V])
plot(sol2, idxs=[test2.src.y])

#=
@mtkmodel PVLoad begin
    @components begin 
        cell = PVCell()
        load = Resistor(R=10)
        i = CurrentSensor()
        g = Ground()
    end

    @equations begin 
        connect(g.g, cell.n, load.n)
        connect(cell.p, i.n)
        connect(i.p, load.p)
        cell.G ~ max(1000*cos(2*pi/10 * t), 0.0)
    end
end

@mtkbuild loaded = PVLoad()
prob = ODEProblem(loaded, [], (0.0, 10.0); guesses=[loaded.load.i => 0])
sol = solve(prob; dtmax=0.01)
=#

@mtkmodel MPPTCell begin
    @parameters begin
        over_v(t) = 0
    end
    
    @components begin 
        cell = PVCell() # determines the voltage draw
        # vref = VoltageSource()
        # res = Resistor(R=2)
        res = BetaMPPTLoad() # determines the current draw given the voltage draw
        i = CurrentSensor()
        # src = BlockComponents.Ramp(;start_time=0, offset=0, height=35, duration=10)
        ground = Ground()
        
        temp = TempSensor()
    end

    @equations begin 
        connect(ground.g, cell.n, res.n)
        connect(cell.p, i.n)
        connect(i.p, res.p)
        connect(cell.Vt, res.Vt)
        # vref.V ~ src.y
        # vref.V ~ t*35/10

        # cell.G ~ max(1000*sin(2*pi*t/10), 0)
        # cell.T_reading ~ 25*sin(2*pi*t/2) + 323

        cell.G ~ temp.G_eff
        cell.T_reading ~ max(temp.T, 125)
    end
    # @continuous_events begin
    #     ModelingToolkit.SymbolicContinuousCallback(;eqs=[cell.V.v ~ cell.Vocn + 0.1], affect=ModelingToolkit.ImperativeAffect(modified=(;over_v, i = cell.Im.I), observed=(;v = cell.V.v, Vocn = cell.Vocn)) do m,o,c,i 
    #         @show o i.t
    #         return (;over_v = o.v >= o.Vocn ? 1.0 : 0.0)
    #     end, affect_neg=Equation[], reinitializealg=OrdinaryDiffEq.BrownFullBasicInit())

    #     ModelingToolkit.SymbolicContinuousCallback(;eqs=[cell.V.v ~ cell.Vocn], affect=Equation[], affect_neg=ModelingToolkit.ImperativeAffect(modified=(;over_v, i = cell.Im.I), observed=(;v = cell.V.v, Vocn = cell.Vocn)) do m,o,c,i 
    #         @show o i.t
    #         return (;over_v = 0)
    #     end, reinitializealg=OrdinaryDiffEq.BrownFullBasicInit())
    # end
end

@mtkbuild cell = MPPTCell()
prob2 = ODEProblem(cell, [cell.res.β => -28], (0.0, end_time); guesses=[cell.cell.Im.i => 8.207600054307171, cell.res.v => 1])
# prob2 = ODEProblem(cell, [], (0.0, end_time); guesses=[cell.cell.Im.i => 8.207600054307171, cell.res.v => 1])
sol2 = solve(prob2; dtmax=0.001)


# plot(sol2[test2.vref.V], sol2[-test2.i.i], title="I-V")
plot(sol2, idxs=[cell.cell.T_reading])
plot(sol2, idxs=[cell.cell.G])
plot(sol2, idxs=[cell.temp.G_eff])
plot(sol2, idxs=[cell.cell.Im.I])
plot(sol2, idxs=[cell.cell.V.v])
plot(sol2, idxs=[cell.cell.Im.I * cell.cell.V.v], title="Power")
plot(sol2.t[1:63], sol2[cell.cell.Im.I][1:63] .* sol2[cell.cell.V.v][1:63])

# Day 1 data
plot(sol2.t[1:1016], sol2[cell.cell.T_reading][1:1016])
plot(sol2.t[1:1016], sol2[cell.cell.G][1:1016])
plot(sol2.t[1:1016], sol2[cell.temp.G_eff][1:1016])
plot(sol2.t[1:1016*29], sol2[cell.cell.Im.I][1:1016*29])


plot(sol2, idxs=[cell.cell.V.v])
plot(sol2, idxs=[cell.vref.V])
plot(sol2, idxs=[cell.src.y])