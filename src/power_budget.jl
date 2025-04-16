using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using Plots

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

facing_sun = Int.(cos.(theta) .>= 0) # 1 if solar panel is facing Sun, 0 otherwise
facing_sun_interp = QuadraticInterpolation(facing_sun, times_adj)


################ 5. Power calculation ################
# Use solar flux and the panel angle to estimate power collected over time.

# @mtkmodel SolarPanel begin
#     @parameters begin
#         G = 1361      # solar irradiance at LEO orbit [W/m^2]
#         A = 5            # solar panel area [m^2]
#         η_ref = 0.3        # reference efficiency
#         T_ref = 300        # reference temperature [K]
#         β = 0.004            # temperature coefficient [1/K]

#         α = 0.9            # absorptivity of panel
#         ϵ = 0.8            # emissivity of panel
#         σ = 5.67e-8       # Stefan-Boltzmann constant [W/(m^2*K^4)]

#     end
#     @variables begin
#         G_eff(t)    # effective irradiance at solar panel [W/m^2]
#         T(t)        # panel temperature [K]
#         η(t)        # panel efficiency
#         P(t)        # power produced [W]
#         θ(t)    # angle between solar panel normal and Sun direction [rad]
#         in_sunlight(t)  # 1 if satellite is in sunlight, 0 otherwise
#         sun_facing(t)   # 1 if solar panel is facing Sun, 0 otherwise

#     end
#     @components begin
#         ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
#         ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
#         ex_facing_sun = TimeVaryingFunction(f=facing_sun_interp)    # 1 if solar panel is facing Sun, 0 otherwise
#         # MPPT 
#         # DC-DC converter (takes in voltage)
#     end
#     @equations begin
#         θ ~ ex_theta.output.u
#         in_sunlight ~ ex_sunlight.output.u
#         sun_facing ~ ex_facing_sun.output.u

#         T ~ ((α * G_eff)/(ϵ * σ))^(1/4)
#         G_eff ~ G * cos(θ) * in_sunlight * sun_facing
#         η ~ η_ref*(1 - β*(T - T_ref))
#         P ~ G_eff * A * η
#     end
# end 

@mtkmodel PanelModel begin
    @components begin
        ex_theta = TimeVaryingFunction(f=theta_interp)                    # angle between solar panel normal and Sun direction [rad]
        ex_sunlight = TimeVaryingFunction(f=sunlight_interp)      # 1 if satellite is in sunlight, 0 otherwise
        ex_facing_sun = TimeVaryingFunction(f=facing_sun_interp)
        mdl = SolarPanel()
    end
    @equations begin
        connect(ex_theta.y, mdl.θ)
        connect(ex_sunlight.y, mdl.in_sunlight)
        connect(ex_facing_sun.y, mdl.sun_facing)
    end
end 

@mtkbuild mdl = PanelModel()
# @named mdl = SolarPanel()
# cmdl = complete(mdl)
prob = ODEProblem(mdl, [], (0, end_time), [])
sol = solve(prob; saveat=times_adj)

plot(sol.t[1:1728000], sol[mdl.P][1:1728000], label="power", xlabel="time", ylabel="power [W]", title="Power Generated", size = (1000,1000))
plot(sol.t[1:1728000], sol[mdl.P][1:1728000], label="power", xlabel="time", ylabel="power [W]", title="Power Generated", size = (1000,1000)) 
plot(sol.t[1:86400], sol[mdl.P][1:86400], label="power", xlabel="time", ylabel="power [W]", title="Power Generated")

plot(sol.t[1:43200], sol[mdl.theta][1:43200], label="θ")
plot(sol.t[1:43200], sol[mdl.sun_facing][1:43200], label="sun_facing")
plot!(sol.t[1:43200], sol[cos(mdl.θ)][1:43200], label="cos(θ)")
plot(sol.t[1:43200], sol[mdl.in_sunlight][1:43200], label="sunlight")
plot(sol.t[1:3200], sol[mdl.G_eff][1:3200], label="G_eff")