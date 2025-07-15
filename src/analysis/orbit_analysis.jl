# This is the orbital analysis used in main.jl
using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations

################ 1. Orbital mechanics ################
# Compute the satellite's position over time using orbital parameters.

jd₀ = date_to_jd(2026, 1, 1, 0, 0, 0)  # start time
days = 30  

orb = KeplerianElements(
           jd₀,                 # epoch [JD]
           7171e3,              # semi-major axis [m]
           0.001,               # small eccentricity
           60   |> deg2rad,     # inclination [rad]
           100    |> deg2rad,   # RAAN [rad]
           90     |> deg2rad,   # argument of perigee [rad]
           20     |> deg2rad    # true anomaly [rad]
       )

orbp = Propagators.init(Val(:J2), orb)

ret = Propagators.propagate!(orbp, 0:1:(86400*days)) # propagate for 30 days

sat_pos = getindex(ret,1) # collect the position data
sat_vel = getindex(ret,2) # collect the velocity data


################ 2. Sun vector ################
# Determine where the sun is in relation to the satellite at any given time.

Δt = 1.0    # 1-minute time step in seconds
times = collect(jd₀:Δt/86400:jd₀ + days)  # 86400 = seconds per day
times_adj = times .- jd₀ # adjusted time array where t=0 corresponds to jd₀ 
end_time = maximum(times_adj)

sun_pos_mod = [sun_position_mod(jd) for jd in times] # vector from Earth center to Sun
sun_pos_teme = @. r_eci_to_eci((MOD(),), times, (TEME(),), times) * sun_pos_mod
sun_vec = sun_pos_teme .- sat_pos # vector from Satellite to Sun


################ 3. Illumination condition ################
# Figure out whether the satellite is in sunlight or eclipse.

sat_condition = lighting_condition.(sat_pos, sun_pos_teme)
sunlight = Int.(sat_condition .== :sunlight)   # 1 if satellite is in sunlight, 0 otherwise  
sunlight_interp =  QuadraticInterpolation(sunlight, times_adj)    


################ 4. Solar panel orientation ################
# Compute the angle between the panel normal and the sun vector. 

sun_unit = sun_vec./norm.(sun_vec) # unit vector from Satellite to Sun
panel_norm = sat_vel./norm.(sat_vel) # unit norm for solar panels; assume panels face direction of velocity

theta = acos.(dot.(panel_norm, sun_unit)) # angle between solar panel normal and Sun direction [rad]
theta_interp =  QuadraticInterpolation(theta, times_adj)