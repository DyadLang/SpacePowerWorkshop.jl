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
        sun_facing(t)   # 1 if solar panel is facing Sun, 0 otherwise
    end
    @components begin
        ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
        ex_facing_sun = TimeVaryingFunction(f=facing_sun_interp)    # 1 if solar panel is facing Sun, 0 otherwise
        # MPPT 
        # DC-DC converter (takes in voltage)
    end
    @equations begin
        θ ~ ex_theta.output.u
        in_sunlight ~ ex_sunlight.output.u
        sun_facing ~ ex_facing_sun.output.u

        T ~ max(((α * G_eff)/(ϵ * σ))^(1/4), 2.7)
        G_eff ~ G * cos(θ) * in_sunlight * sun_facing
    end
end 

@mtkbuild temp = TempSensor() 
prob = ODEProblem(temp, [], (0, end_time), [])
sol = solve(prob;  saveat=times_adj)

plot(sol.t[1:43200], sol[temp.T][1:43200], label="temp.T")
plot(sol.t[1:43200], sol[temp.G_eff][1:43200], label="temp.G_eff")
plot(sol.t[1:43200], sol[temp.θ][1:43200], label="temp.θ")
plot(sol.t[1:43200], sol[temp.in_sunlight][1:43200], label="temp.in_sunlight")
plot(sol.t[1:43200], sol[temp.sun_facing][1:43200], label="temp.sun_facing")

@mtkmodel PVTest begin
    @components begin 
        cell = PVCell()
        vref = VoltageSource()
        i = CurrentSensor()
        src = BlockComponents.Ramp(;start_time=0, offset=0, height=35, duration=10)
        ground = Ground()

        temp = TempSensor()
    end

    @equations begin 
        connect(ground.g, cell.n, vref.n)
        connect(cell.p, i.n)
        connect(i.p, vref.p)
        vref.V ~ src.y
        # cell.G ~ 1000

        cell.G ~ temp.G_eff
        temp.T ~ cell.T_reading
    end
end

@mtkbuild test = PVTest()
prob = ODEProblem(test, [], (0.0, end_time); guesses=[test.cell.Rs_c.i => 0])
sol = solve(prob)

plot(sol[test.vref.V], sol[test.i.i], title="I-V")

plot(sol, idxs=[test.cell.V.v])
plot!(sol, idxs=[test.vref.V])
plot(sol, idxs=[test.cell.I.i])
plot!(sol, idxs=[test.i.i])
plot(sol, idxs=[test.cell.T_reading])
plot(sol, idxs=[test.cell.G])
plot(sol, idxs=[test.temp.in_sunlight])
plot(sol, idxs=[test.temp.sun_facing])

plot(sol, idxs=[test.temp.T])
plot(sol, idxs=[test.temp.G_eff])
plot(sol, idxs=[test.temp.θ])

plot(sol, idxs=[test.vref.V])
plot(sol, idxs=[test.cell.rolloff])
plot(sol, idxs=[test.cell.ipv])
plot(sol, idxs=[test.cell.Im.I])

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