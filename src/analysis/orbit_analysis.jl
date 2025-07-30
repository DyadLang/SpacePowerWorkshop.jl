# This is the orbital analysis used in main.jl
using SatelliteToolbox
using SatelliteAnalysis
using LinearAlgebra
using DataInterpolations
import ModelingToolkit as MTK

using DyadInterface

export OrbitalTransientAnalysis, OrbitalTransientAnalysisSpec, OrbitalTransientAnalysisSolution, AbstractOrbitalAnalysisSpec

abstract type AbstractOrbitalTransientAnalysisSpec <: DyadInterface.AbstractAnalysisSpec end

Base.@kwdef mutable struct OrbitalTransientAnalysisSpec{M} <: AbstractOrbitalTransientAnalysisSpec
    name::Symbol = :OrbitalTransientAnalysis
    alg::String = "Rodas5P"
    abstol::Real = 1e-6
    reltol::Real = 1e-6
    saveat::Real = 0
    dtmax::Real = 1e-3
    IfLifting::Bool = false

    sunlight_interp_param::String
    theta_interp_param::String

    # The start date of the analysis, **as a Julian date**.
    # Obtain this from `SatelliteToolbox.date_to_jd(year, month, day, hh, mm, ss)`
    start_date::Real = date_to_jd(2026, 1, 1, 0, 0, 0)
    # The number of days to simulate after start_date
    days::Real = 30.0
    # The semi-major axis of the orbit
    semimajor_axis::Real = 7171e3
    # The eccentricity of the orbit (small e)
    eccentricity::Real = 0.001
    # The inclination of the orbit, in degrees
    inclination::Real = 60
    # The RAAN of the orbit, in degrees
    raan::Real = 100
    # The argument of perigee, in degrees
    argument_of_perigee::Real = 90
    # The true anomaly, in degrees
    anomaly::Real = 20
    # The solar panel model to simulate.  This must extend or 
    # satisfy the [`PartialSolarPanel`](@ref) interface.
    model::M
end

struct OrbitalTransientAnalysisSolution{Spec, Sol, T, Ts, Pos, Vel, SunPos, SI <: DataInterpolations.AbstractInterpolation, TI <: DataInterpolations.AbstractInterpolation} <: DyadInterface.AbstractAnalysisSolution
    spec::Spec
    "The solution to the ODE problem."
    sol::Sol
    "The start date of the analysis as a Julian date.  Add this to all elements of `time_rel` to get the absolute time in Julian days."
    start_date::T
    "The time relative to the start date as Julian days."
    time_rel::Ts
    "The satellite position in the TEME reference frame."
    sat_pos::Pos
    "The satellite velocity in the TEME reference frame."
    sat_vel::Vel
    "The sun position in the TEME reference frame."
    sun_pos::SunPos
    "The sunlight condition interpolation (callable with time as a parameter)."
    sunlight::SI
    "The solar panel orientation interpolation (callable with time as a parameter)."
    theta::TI
end

function Base.show(io::IO, ::MIME"text/plain", sol::OrbitalTransientAnalysisSolution)
    printstyled(io, "Orbital Transient Analysis Solution", color=:cyan, bold=true)
    println(io, " from $(julian2datetime(sol.spec.start_date)) to $(julian2datetime(sol.start_date + sol.time_rel[end]))")
    println(io, " $(length(sol.time_rel)) time steps")
    println(io, " with fields time_rel, sat_pos, sat_vel, sun_pos, sunlight, theta")
end

OrbitalTransientAnalysis(; kwargs...) = DyadInterface.run_analysis(OrbitalTransientAnalysisSpec(; kwargs...))

# TODO better interface for this, why should I have to fake this out with a transient?!
function DyadInterface.ODEProblemConfig(spec::OrbitalTransientAnalysisSpec)
    transient_spec = DyadInterface.TransientAnalysisSpec(;
        alg = spec.alg,
        abstol = spec.abstol,
        reltol = spec.reltol,
        saveat = spec.saveat,
        dtmax = spec.dtmax,
        IfLifting = spec.IfLifting,
        start = 0.0, 
        stop = spec.days,
        model = nothing
    )
    return DyadInterface.ODEProblemConfig(transient_spec)
end

function DyadInterface.run_analysis(spec::OrbitalTransientAnalysisSpec)
    # ## 1. Orbital mechanics.
    # Compute the satellite's position over time using orbital parameters.
    # See directly below this for the implementation of this function.
    (; time_rel, sat_pos, sat_vel, sun_pos_teme, sunlight_interp, theta_interp) = run_orbit(;
        start_date = spec.start_date,
        days = spec.days,
        semimajor_axis = spec.semimajor_axis,
        eccentricity = spec.eccentricity,
        inclination = spec.inclination,
        raan = spec.raan,
        argument_of_perigee = spec.argument_of_perigee,
        anomaly = spec.anomaly,
    )
    # ## 2. Set up the model.
    config = DyadInterface.ODEProblemConfig(spec)
    simplified_model = DyadInterface.get_simplified_model(spec)
    prob = DyadInterface.setup_prob(
        simplified_model, 
        [
            MTK.parse_variable(simplified_model, spec.theta_interp_param) => theta_interp, 
            MTK.parse_variable(simplified_model, spec.sunlight_interp_param) => sunlight_interp,
            simplified_model.converter.β => -28, 
            simplified_model.converter.stored_energy => 1,
        ], 
        config.tspan, 
        (;
            guesses = [simplified_model.cell.Im.i => 8.207600054307171, simplified_model.converter.v => 1]
        )
    )

    sol = solve(prob, config.alg; config.saveat, config.abstol, config.reltol, config.dtmax)

    return OrbitalTransientAnalysisSolution(
        spec, 
        sol,
        spec.start_date,
        time_rel,
        sat_pos,
        sat_vel,
        sun_pos_teme,
        sunlight_interp,
        theta_interp,
    )
end

# function DyadInterface.run_analysis(spec::OrbitalGUI)
#     res = run_analysis(OrbitalTransientAnalysisSpec(; spec...))
#     gui = build_figure(res)
#     screen = display(gui)
#     wait(screen)
#     return res
# end

function run_orbit(; start_date, days, semimajor_axis, eccentricity, inclination, raan, argument_of_perigee, anomaly)

    orb = KeplerianElements(
            start_date,                       # epoch [JD]
            semimajor_axis,                   # semi-major axis [m]
            eccentricity,                     # small eccentricity
            inclination         |> deg2rad,   # inclination [rad]
            raan                |> deg2rad,   # RAAN [rad]
            argument_of_perigee |> deg2rad,   # argument of perigee [rad]
            anomaly             |> deg2rad    # true anomaly [rad]
        )

    orbp = Propagators.init(Val(:J2), orb)

    ret = Propagators.propagate!(orbp, 0:1:(86400*days)) # propagate for 30 days

    # Note: these are both in the TEME frame....
    sat_pos = getindex(ret,1) # collect the position data
    sat_vel = getindex(ret,2) # collect the velocity data


    ################ 2. Sun vector ################
    # Determine where the sun is in relation to the satellite at any given time.

    Δt = 1.0    # 1-second time step in seconds
    times = collect((start_date):Δt/86400:(start_date + days))  # 86400 = seconds per day
    times_adj = times .- start_date # adjusted time array where t=0 corresponds to jd₀ 
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

    return (; time_rel = times_adj, sat_pos, sat_vel, sun_pos_teme, sunlight_interp, theta_interp)
end

(; time_rel, sat_pos, sat_vel, sun_pos_teme, sunlight_interp, theta_interp) = run_orbit(;
    start_date = date_to_jd(2026, 1, 1, 0, 0, 0),
    days = 30,
    semimajor_axis = 7171e3,
    eccentricity = 0.001,
    inclination = 60,
    raan = 100,
    argument_of_perigee = 90,
    anomaly = 20,
)