using GLMakie, GeoMakie

using SatelliteToolbox, SatelliteToolboxTransformations
using DataInterpolations: LinearInterpolation

using SpacePowerWorkshop
using SpacePowerWorkshop: jd_to_gmst, date_to_jd

using LinearAlgebra, Dates

# ## Setup
# We will use some image assets for the visualization, these should all be publically usable.

# Full-sky mosaic from the European Southern Observatory.
# [Source webpage here.](https://www.eso.org/public/images/eso0932a/)
if !isfile("eso0932a.tif")
    download("https://cdn.eso.org/images/original/eso0932a.tif", "eso0932a.tif")
end
# NASA full-earth blue marble in July.
# [acknowledgment](https://visibleearth.nasa.gov/images/74092/july-blue-marble-next-generation)
if !isfile("bluemarble.png")
    download("https://eoimages.gsfc.nasa.gov/images/imagerecords/76000/76487/world.200406.3x5400x2700.png", "bluemarble.png")
end
# Load these images into memory as matrices of RGB
using FileIO, ImageIO
blue_marble_img = load("bluemarble.png")
starry_background_img = load("eso0932a.tif")

# Finally, we run a simulation of the solar cell so we can show a dashboard of the power output.
using ModelingToolkit, DyadInterface
@time @mtkbuild panel = SolarPanel()
prob = ODEProblem(
    panel, 
    [panel.converter.β => -28, panel.converter.stored_energy => 1], # parameter overrides
    (0.0, SpacePowerWorkshop.end_time);                             # time span
    guesses = [panel.cell.Im.i => 8.207600054307171, panel.converter.v => 1]
)
sol = solve(prob; dtmax=0.001)

# ### Satellite orbit
# We want everything to be in the reference frame ECEF (Earth-Centered, Earth-Fixed)
# since that is the best for GIS visualization.
# We could choose to keep everything in TEME which is inertial and so the Earth rotates within it,
# but that gets annoying to deal with until we have better tools for it in GeoMakie (which we will at some point!)
#
# So, let's transform everything from TEME to ECEF.
# The propagator is in TEME.
eop = SatelliteToolboxTransformations.fetch_iers_eop()
sat_pos_ecef = @. SatelliteToolboxTransformations.r_eci_to_ecef(
    (TEME(),),
    (ITRF(),),
    SpacePowerWorkshop.times, # times in julian date
    (eop,)
) * SpacePowerWorkshop.sat_pos

# sun_light = Makie.DirectionalLight(Colors.WP_A |> RGBf, direction_from_sun_to_earth)
fig = with_theme(theme_dark()) do
    Figure(; figure_padding = 0, size = (600, 900));
end
ax = GlobeAxis(fig[1, 1]; show_axis = false, height = 600, width = 600)
# push!(ax.scene.lights, sun_light)
earth_plt = meshimage!(
    ax, -180..180, -90..90, blue_marble_img;
    uv_transform = :rotr90,
    shading = Makie.MultiLightShading
)
background_plt = meshimage!(
    ax, -180..180, -90..90, starry_background_img; 
    uv_transform = :rotr90, zlevel = 2e8,
    shading = NoShading
)
# BUG SECTION
# 1. The geomakie meshimage recipe does not propagate shading
only(earth_plt.plots).shading[] = Makie.MultiLightShading
# END BUG SECTION

fig

time_rel = Observable(0.001)

satellite_marker = lift(time_rel) do t
    sat_pos_ecef[round(Int, t * 86400) + 1]
end

satellite_trajectory = lift(time_rel) do t
    sat_pos_ecef[max(1, round(Int, t * 86400) - 180 * 30):round(Int, t * 86400) + 1]
end

sunlight_fraction = lift(time_rel) do t
    SpacePowerWorkshop.sunlight_interp(t)
end

satellite_marker_plt = scatter!(
    ax.scene,
    @lift([$satellite_marker]);
    # marker = satellite_mesh,
)

satellite_trajectory_spec = lift(satellite_trajectory) do trajectory
    [Makie.SpecApi.Lines(
        trajectory;
        color = if length(trajectory) == 1
            [Makie.wong_colors(1.0)[2]]
            else
                RGBAf.((Makie.wong_colors()[2],), LinRange(0, 1, length(trajectory)))
            end
    )]
end
satellite_trajectory_plt = plotlist!(ax.scene, satellite_trajectory_spec)



# ## Animation
# Since we set up this observable update chain, 
# updating the time observable will cause a cascade of updates
# down to the rest of the observables.
@time record(fig, "earth_rotations.mp4", LinRange(2, 4, 2000); framerate = 60, update = false) do t_days
    t_seconds = t_days * 86400
    # Update the relative time Observable
    # This cascades and updates all of the other observables,
    # which update the visuals.
    time_rel[] = t_days
    # reset_limits!(power_ax)
end

screen = display(fig; update = false)
wait(screen) # wait for display to close
