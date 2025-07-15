using GLMakie, GeoMakie
using DataInterpolations: LinearInterpolation
using SpacePowerWorkshop
using SpacePowerWorkshop: jd_to_gmst, date_to_jd
using LinearAlgebra

using SatelliteToolboxTransformations

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


# satellite_mesh = load("/Users/anshul/Downloads/ImageToStl.com_CubeSat+-+1+RU+Generic/CubeSat - 1 RU Generic.obj")

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

# Finally, create a 'dashboard layout' at the bottom of the figure,'
# that we can place information in.
info_gl = GridLayout(fig[1, 1, Bottom()], height = 100, tellheight = false, valign = :bottom, alignmode = Outside(10))

ind_grad = Makie.cgrad(:RdYlGn_10)
in_sunlight_label = Label(
    info_gl[1, 1], 
    @lift Makie.rich("Sunlight fraction: ", Makie.rich("⬤", color = ind_grad[$sunlight_fraction]), string(round($sunlight_fraction, digits = 2)));
    valign = :center,
    tellheight = false,
)

min_power, max_power = extrema(sol[panel.cell.Im.I * panel.cell.V.v])

power_label = Label(
    info_gl[1, 2],
    @lift Makie.rich("Power: ", 
        Makie.rich(
            lpad(string(round(Int, (sol($(time_rel); idxs = panel.cell.Im.I * panel.cell.V.v) - min_power) / (max_power - min_power) * 100)), 3);
            color = Makie.wong_colors(1.0)[1],
            font = "Fira Mono"
            ), 
        "%",
    );
    valign = :center,
    tellheight = false,
)
power_ax, power_plot = lines(
    info_gl[1, 3],
    lift(time_rel) do t
        trange = LinRange(max(0, t - 90/(60*24)), t, 1000)
        current_power = sol(trange; idxs = panel.cell.Im.I * panel.cell.V.v)
        Point2f.(LinRange(0, 1, length(trange)), current_power)
    end;
    axis = (; 
        limits = ((0, 1), (min_power, max_power)),
        xticks = ([0, 0.5, 1], ["-90 min", "-45 min", "now"]),
    ),
    linewidth = 5
)

# ## Ground track
# We can also plot a ground track of the satellite.
# Maybe this should have a Tyler map?  But that's too brittle.
# This gets plotted on a projected map of earth.
ground_ax = GeoAxis(fig[2, 1]; limits = ((-180, 180), (-90, 90)), title = "Ground Track")
meshimage!(ground_ax, -180..180, -90..90, blue_marble_img; uv_transform = :rotr90)
lines!(ground_ax, GeoMakie.coastlines(); color = :white)
satellite_marker_ground_plt = scatter!(ground_ax, @lift([$satellite_marker]); source = "+proj=cart +type=crs", marker = :circle, color = :blue, strokecolor = :white, strokewidth = 1)

# this janky thing is because Makie doesn't allow the transformation to 
# manipulate geometry directly.
# but that's coming soon!
trajectory_longlat = lift(satellite_trajectory) do trajectory_ecef
    trajectory_longlat = [begin; lat, lon, alt = SatelliteToolboxTransformations.ecef_to_geodetic(ecef); Point2(rad2deg(lon), rad2deg(lat)); end for ecef in trajectory_ecef]
    split_mls = Base.split(GeometryBasics.LineString(trajectory_longlat), 180.0)
    joined_ls = Point2{Float64}[] 
    # joined_ls_color = Float64[]
    i = 1
    for idx in 1:(length(split_mls.linestrings) - 1)
        linestring = split_mls.linestrings[idx]
        append!(joined_ls, getproperty(linestring, :points))
        # append!(joined_ls_color, collect(Float64(i):Float64(i + length(getproperty(linestring, :points)) - 1)))
        push!(joined_ls, Point2{Float64}(NaN))
        # push!(joined_ls_color, NaN)
        i += length(getproperty(linestring, :points))
    end
    linestring = split_mls.linestrings[end]
    append!(joined_ls, getproperty(linestring, :points))
    # append!(joined_ls_color, collect(Float64(i):Float64(i + length(getproperty(linestring, :points)) - 1)))
    # [Makie.SpecApi.Lines(joined_ls, color = joined_ls_color)]
    joined_ls
end
# finally, plot the ground track.
satellite_split_trajectory_plt = lines!(ground_ax, trajectory_longlat; color = Makie.Cycled(2))


# ## Animation
# Since we set up this observable update chain, 
# updating the time observable will cause a cascade of updates
# down to the rest of the observables.
@time record(fig, "earth_rotations.mp4", LinRange(0.01, 2, 2000); framerate = 60, update = false) do t_days
    t_seconds = t_days * 86400
    # Update the relative time Observable
    # This cascades and updates all of the other observables,
    # which update the visuals.
    time_rel[] = t_days
    # reset_limits!(power_ax)
end


tmax = 1.0
indmax = 86400 # 1 day in seconds
