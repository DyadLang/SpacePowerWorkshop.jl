using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using Plots
using Example1

################ 1. Orbital mechanics ################
# Compute the satellite's position over time using orbital parameters.

jd₀ = date_to_jd(2025, 1, 1, 0, 0, 0)  # start time
days = 30  

orb = KeplerianElements(
           jd₀,                 # epoch [JD]
           7190.982e3,          # semi-major axis [m]
           0.001111,            # small eccentricity
           98.405 |> deg2rad,   # inclination [rad]
           100    |> deg2rad,   # RAAN [rad]
           90     |> deg2rad,   # argument of perigee [rad]
           19     |> deg2rad    # true anomaly [rad]
       )

orbp = Propagators.init(Val(:J2), orb)

ret = Propagators.propagate!.(orbp, 0:1:(86400*days)) # propagate for 30 days

sat_pos = getindex.(ret,1) # collect the position data
sat_vel = getindex.(ret,2) # collect the velocity data


################ 2. Sun vector ################
# Determine where the sun is in relation to the satellite at any given time.

Δt = 1.0    # 1-minute time step in seconds
times = collect(jd₀:Δt/86400:jd₀ + days)  # 86400 = seconds per day
times_adj = times .- jd₀ # adjusted time array where t=0 corresponds to jd₀ 
end_time = maximum(times_adj)

sun_pos = [sun_position_mod(jd) for jd in times] # vector from Earth center to Sun
sun_vec = sun_pos .- sat_pos # vector from Satellite to Sun

################ 3. Illumination condition ################
# Figure out whether the satellite is in sunlight or eclipse.

sat_condition = lighting_condition.(sat_pos, sun_pos)
sunlight = Int.(sat_condition .== :sunlight)   # 1 if satellite is in sunlight, 0 otherwise  
sunlight_interp =  QuadraticInterpolation(sunlight, times_adj)    

################ 4. Solar panel orientation ################
# Compute the angle between the panel normal and the sun vector. 

sun_unit = sun_vec./norm.(sun_vec) # unit vector from Satellite to Sun
panel_norm = sat_vel./norm.(sat_vel) # unit norm for solar panels; assume panels face direction of velocity

theta = acos.(dot.(panel_norm, sun_unit)) # angle between solar panel normal and Sun direction [rad]
theta_interp =  QuadraticInterpolation(theta, times_adj)


################ 5. Power calculation ################
# Use solar flux and the panel angle to estimate power collected over time.

@mtkmodel PanelModel begin
    @components begin
        ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
        solar_panel = SolarPanelSimple()
    end
    @equations begin
        connect(ex_theta.output.u, solar_panel.θ)
        connect(ex_sunlight.output.u, solar_panel.in_sunlight)
    end
end 

@mtkbuild mdl = PanelModel() 
prob = ODEProblem(mdl, [], (0, end_time), [])
sol = solve(prob; saveat=times_adj)

period = 1 * 86400

plot(sol.t[1:period], sol[mdl.solar_panel.P][1:period], label=false, xlabel="Days", ylabel="Power [W]", title="Power Generated", linewidth=2, size=(800,600))
plot(sol.t[1:period], sol[mdl.solar_panel.T][1:period], label="panel.T")
plot(sol.t[1:period], sol[mdl.solar_panel.G_eff][1:period], label=false, xlabel="Days", ylabel="Irradiance [W/m2]", title="Effective Irradiance", linewidth=2,size=(800,600))
plot(sol.t[1:period], sol[mdl.solar_panel.η][1:period], label="panel.η")
plot(sol.t[1:period], sol[mdl.solar_panel.θ][1:period], label="panel.θ")












################ 6. Old stuff ################

#=
@mtkmodel SolarPanel begin
    @parameters begin
        G = 1361      # solar irradiance at LEO orbit [W/m^2]
        A = 5            # solar panel area [m^2]
        η_ref = 0.3        # reference efficiency
        T_ref = 300        # reference temperature [K]
        β = 0.004            # temperature coefficient [1/K]

        α = 0.9            # absorptivity of panel
        ϵ = 0.8            # emissivity of panel
        σ = 5.67e-8       # Stefan-Boltzmann constant [W/(m^2*K^4)]

    end
    @variables begin
        G_eff(t)    # effective irradiance at solar panel [W/m^2]
        T(t)        # panel temperature [K]
        η(t)        # panel efficiency
        P(t)        # power produced [W]
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

        T ~ ((α * G_eff)/(ϵ * σ))^(1/4)
        G_eff ~ G * cos(θ) * in_sunlight * sun_facing
        η ~ η_ref*(1 - β*(T - T_ref))
        P ~ G_eff * A * η
    end
end 
=#

@mtkbuild solar_panel = SolarPanel() 
# prob = ODEProblem(mdl, [], (0, end_time), []; guesses=[mdl.resistor.i => 0])
prob = ODEProblem(solar_panel, [], (0, end_time), [])
# prob = ODEProblem(mdl, [], (0, end_time), []; guesses=[mdl.panel.v => 15, mdl.panel.Ir => 0, mdl.resistor.i => 0])
sol = solve(prob;  saveat=times_adj)
plot(sol.t[1:43200], sol[solar_panel.θ][1:43200])


@mtkmodel PVArrayVerification begin
    @components begin
        Gn = BlockComponents.Const(k=1000)
        Tn = BlockComponents.Const(k=298.15)
        source = VoltageSource()
        ground = Ground()
        pVArray = PVArray()
        rampVoltage = BlockComponents.Ramp(duration=1, height=45, offset=-10)

        ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
        ex_facing_sun = TimeVaryingFunction(f=facing_sun_interp)
    end
    @equations begin
        connect(rampVoltage.y, source.V)
        connect(Gn.y, pVArray.G)
        connect(Tn.y, pVArray.T)
        connect(pVArray.p, source.p)
        connect(pVArray.n, source.n)
        connect(ground.g, source.n)

        connect(ex_theta.output.u, pVArray.θ)
        connect(ex_sunlight.output.u, pVArray.in_sunlight)
        connect(ex_facing_sun.output.u, pVArray.sun_facing)
    end

end

@mtkmodel PanelModel begin
    @components begin
        ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
        ex_facing_sun = TimeVaryingFunction(f=facing_sun_interp)
        # panel = SolarPanel()
        # panel = CurrentSource()
        panel = PVArray()
        ground = Ground()
        resistor = Resistor(R=10)
    end
    @equations begin
        connect(ex_theta.output.u, panel.θ)
        connect(ex_sunlight.output.u, panel.in_sunlight)
        connect(ex_facing_sun.output.u, panel.sun_facing)
        connect(panel.n, ground.g, resistor.n)
        connect(panel.p, resistor.p)
    end
end 

@mtkbuild mdl = PanelModel() 
# prob = ODEProblem(mdl, [], (0, end_time), []; guesses=[mdl.resistor.i => 0])
prob = ODEProblem(mdl, [], (0, end_time), []; guesses=[mdl.panel.v => 15, mdl.panel.Ir => 0, mdl.panel.Ipv => 1, mdl.resistor.i => 1])
# prob = ODEProblem(mdl, [], (0, end_time), []; guesses=[mdl.panel.v => 15, mdl.panel.Ir => 0, mdl.resistor.i => 0])
sol = solve(prob; dtmin = 0.00001, saveat=times_adj)

plot(sol.t[1:1728000], sol[mdl.panel.P][1:1728000], label="power", xlabel="time", ylabel="power [W]", title="Power Generated", size = (1000,1000))
plot(sol.t[1:86400], sol[mdl.panel.P][1:86400], label="power", xlabel="time", ylabel="power [W]", title="Power Generated")

plot(sol.t[1:86400], sol[mdl.panel.i * mdl.panel.v][1:86400], label="power", xlabel="time", ylabel="power [W]", title="Power Generated")
plot(sol.t[1:40000], sol[mdl.panel.T][1:40000], label="panel.T")
plot!(sol.t[1:43200], sol[mdl.panel.G/200][1:43200], label="panel.G")
plot(sol.t, sol[mdl.panel.P], label="panel.P")
plot(sol.t[1:43200], sol[mdl.panel.v][1:43200], label="panel.v")
plot(sol.t[1:43200], sol[mdl.panel.i][1:43200], label="panel.i")
plot(sol.t[1:43200], sol[(mdl.panel.v - mdl.panel.Rs*mdl.panel.i)/mdl.panel.a/50][1:43200], label="panel.")
plot(sol.t, sol[mdl.panel.i_vneg], label="panel.i_vneg")
plot!(sol.t, sol[mdl.panel.i_vpos], label="panel.i_vpos")
plot!(sol.t[1:43200], sol[mdl.panel.Ipv][1:43200], label="panel.Ipv")
plot(sol.t[1:40000], sol[mdl.panel.I0][1:40000], label="panel.I0")
plot(sol.t[1:40000], sol[mdl.panel.Id][1:40000], label="panel.Id")
plot(sol.t[1:40000], sol[mdl.panel.Ir][1:40000], label="panel.Ir")
plot(sol.t[1:40000], sol[mdl.panel.Vt][1:40000], label="panel.Vt")
plot!(sol.t[1:40000], sol[mdl.resistor.i][1:40000], label="resistor.i")
plot(sol.t[1:43200], sol[mdl.resistor.v][1:43200], label="resistor.v")

plot(sol.t, sol[mdl.panel.Ir], label="panel.Ir")
plot(sol.t[1:3600], sol[mdl.panel.i][1:3600], label="panel.i")
plot(sol.t, sol[mdl.resistor.i], label="resistor.i")
plot(sol.t[1:3600], sol[mdl.panel.Iph][1:3600], label="Iph")
plot(sol.t[1:3600], sol[mdl.panel.Is][1:3600], label="Is")
plot!(sol.t[1:3600], sol[mdl.panel.ex_term][1:3600], label="ex_term")
plot(sol.t[1:3600], sol[mdl.panel.Iph / ((mdl.panel.ex_term-1) * mdl.panel.Is)][1:3600], label="calculated i")

plot(sol.t[1:3600], sol[mdl.panel.Iph][1:3600], label="Iph")
plot!(sol.t[1:3600], sol[(mdl.panel.ex_term-1) * mdl.panel.Is][1:3600], label="other part")


plot!(sol.t, sol[mdl.panel.sun_facing], label="sun_facing")
plot!(sol.t, sol[cos(mdl.panel.θ)], label="cos(θ)")
plot(sol.t[1:3600], sol[mdl.panel.G_eff][1:3600], label="G_eff")

plot(sol.t[1:43200], sol[mdl.panel.θ][1:43200], label="θ")
plot(sol.t[1:43200], sol[mdl.panel.sun_facing][1:43200], label="sun_facing")
plot!(sol.t[1:43200], sol[cos(mdl.panel.θ)][1:43200], label="cos(θ)")
plot(sol.t[1:43200], sol[mdl.panel.in_sunlight][1:43200], label="sunlight")
plot(sol.t[1:3200], sol[mdl.panel.G_eff][1:3200], label="G_eff")

@mtkmodel TestThing begin
    @components begin
        source = CurrentSource()
        ground = Ground()
        resistor = Resistor(R=10)
    end
    @equations begin
        connect(source.n, ground.g, resistor.n)
        connect(source.p, resistor.p)
    end
end

@mtkbuild thing = TestThing()
prob2 = ODEProblem(thing, [], (0, end_time), []; guesses=[mdl.resistor.i => 0])
sol2 = solve(prob2; saveat=times_adj)