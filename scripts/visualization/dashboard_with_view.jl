using GLMakie
using GeoMakie
using GeoMakie.Proj, GeoMakie.GeometryBasics
using GeoMakie.Makie.Observables

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
prob = @time ODEProblem(
    panel, 
    [panel.converter.β => -28, panel.converter.stored_energy => 1], # parameter overrides
    (0.0, SpacePowerWorkshop.end_time);                             # time span
    guesses = [panel.cell.Im.i => 8.207600054307171, panel.converter.v => 1]
)
sol = @time solve(prob; dtmax=0.001)

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
    Figure(; figure_padding = 0, size = (1000, 600));
end
ax = GlobeAxis(fig[1, 1]; show_axis = false, width = 600, halign = :left, tellheight = true, tellwidth = true)
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

# ### View point for overview plot
# We want to look down at the earth from some point in space,
# usually centered on a point on the earth.  Since JuliaCon is coming up,
# let's look at Pittsburgh.
pittsburgh_lonlat = (-79.9428, 40.4432)
pittsburgh_ecef = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.wgs84)(
        GeoMakie.Geodesy.LLA(; 
        lon = pittsburgh_lonlat[1], 
        lat = pittsburgh_lonlat[2], 
        alt = 2e7
    )
)
# Now, we update the camera to look at Pittsburgh.
cc = cameracontrols(ax.scene)
cc.eyeposition[] = pittsburgh_ecef
cc.lookat[] = Vec3d(0,0,0)
cc.upvector[] = Vec3d(0,0,1)
Makie.update_cam!(ax.scene, cc)

# ## Satellite animation and Observables
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

sun_line = lift(time_rel) do t
    idx = round(Int, t * 86400) + 1
    sat_pos = sat_pos_ecef[idx]
    sun_pos = r_eci_to_ecef(TEME(), ITRF(), t + SpacePowerWorkshop.jd₀, eop) * SpacePowerWorkshop.sun_pos_teme[idx]
    Point3d[sun_pos, sat_pos]
end

sun_line_plt = lines!(ax, sun_line; source = "+proj=cart +type=crs", color = :yellow, linewidth = 2, xautolimits = false, yautolimits = false, reset_limits = false)
on(sunlight_fraction) do sunlit
    sun_line_plt.visible[] = sunlit > 0.1
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
info_gl = GridLayout(
    fig[1, 1, Bottom()], 
    height = 100, 
    tellheight = false, 
    valign = :bottom, 
    alignmode = Outside(10)
)

ind_grad = Makie.cgrad(:RdYlGn_10)
in_sunlight_label = Label(
    info_gl[1, 1], 
    @lift Makie.rich("Sunlight fraction: ", Makie.rich("⬤", color = ind_grad[$sunlight_fraction]), string(round($sunlight_fraction, digits = 2)));
    valign = :center,
    tellheight = true,
)

min_battery, max_battery = extrema(sol[panel.converter.stored_energy])

battery_label = Label(
    info_gl[1, 2],
    @lift Makie.rich("Battery: ", 
        Makie.rich(
            lpad(string(round(Int, (sol($(time_rel); idxs = panel.converter.stored_energy) - min_battery) / (max_battery - min_battery) * 100)), 3);
            color = Makie.wong_colors(1.0)[1],
            font = "Fira Mono"
            ), 
        "%",
    );
    valign = :center,
    tellheight = true,
)

time_label = Label(
    info_gl[2, 1:2],
    lift(time_rel) do t
        Makie.rich(
            "Time: ",
            Dates.format(julian2datetime(t + SpacePowerWorkshop.jd₀), "u dd HH:MM:SS")
        )
    end;
    tellwidth = false,
    tellheight = true,
)

battery_ax, battery_plot = lines(
    info_gl[1:2, 3],
    lift(time_rel) do t
        trange = LinRange(max(0, t - 90/(60*24)), t, 1000)
        current_battery = sol(trange; idxs = panel.converter.stored_energy)
        Point2f.(LinRange(0, 1, length(trange)), current_battery)
    end;
    axis = (; 
        limits = ((0, 1), (min_battery, max_battery)),
        xticks = ([0, 0.5, 1], ["-90 min", "-45 min", "now"]),
        alignmode = Outside(), # important so that the protrusions remain within the grid cell
    ),
    linewidth = 5
)

# ## Ground track
# We can also plot a ground track of the satellite.
# Maybe this should have a Tyler map?  But that's too brittle.
# This gets plotted on a projected map of earth.
diag_gl = GridLayout(fig[1, 2]; alignmode = Outside())
ground_ax = GeoAxis(diag_gl[1, 1]; limits = ((-180, 180), (-90, 90)), title = "Ground Track")
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

satellite_position_longlatalt = lift(satellite_marker) do ecef
    lat, lon, alt = SatelliteToolboxTransformations.ecef_to_geodetic(ecef)
    Point3d(rad2deg(lon), rad2deg(lat), alt)
end

projection = lift(satellite_position_longlatalt) do pos
    "+proj=geos +h=$(pos[3]) +lon_0=$(pos[1]) +lat_0=$(pos[2]) +type=crs"
end

view_label = Label(diag_gl[2, 1], "Satellite View"; halign = :center, font = :bold,tellheight = true, tellwidth = false)
view_ax = GlobeAxis(diag_gl[3, 1]; show_axis = false)
meshimage!(view_ax, -180..180, -90..90, blue_marble_img; uv_transform = :rotr90)

on(satellite_marker) do ecef
    cc = Makie.cameracontrols(view_ax.scene)
    cc.lookat[] = Vec3d(0,0,0)
    cc.eyeposition[] = ecef
    Makie.update_cam!(view_ax.scene, cc)
end



# Here's the dashboard controls to run the animation interactively.
play_button = Button(diag_gl[4, 1]; tellwidth = false, tellheight = true, label = "▶")
is_playing = Observable(false)

play_button_listener = on(play_button.clicks; priority = 1000) do _
    is_playing[] = !is_playing[]
    if is_playing[]
        play_button.label[] = "||"
    else
        play_button.label[] = "▶"
    end
end

player_listener = Makie.Observables.on(events(fig).tick) do tick
    if is_playing[]
        tic = time()
        time_rel[] += 1/3000
        yield()
        toc = time()
    else
        # do nothing
    end
end




@time record(fig, "dashboard_with_view.mp4", LinRange(3, 3.5, 3000); framerate = 60, update = false) do t
    time_rel[] = t
end