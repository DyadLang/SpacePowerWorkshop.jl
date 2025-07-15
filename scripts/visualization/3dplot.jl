using GLMakie, GeoMakie
using DataInterpolations: LinearInterpolation
using SpacePowerWorkshop
using SpacePowerWorkshop: jd_to_gmst, date_to_jd
using LinearAlgebra

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

# ### Sun light
# Since this demo is about solar power, we want to 
# indicate where the sun is with respect to the Earth.
# In order to do this, we use a directional light to indicate the direction of the sun.

# To do this, we need to know the direction of the sun from the Earth.
# We can get this from the sun vector in the propagator.
direction_from_earth_to_sun = normalize(first(SpacePowerWorkshop.sun_vec))
direction_from_sun_to_earth = -direction_from_earth_to_sun

# Rotate the sun direction by 23 degrees south to position it on the equator
# This accounts for Earth's axial tilt
obliquity_angle = deg2rad(23.0)  # 23 degrees in radians
# Create rotation axis (perpendicular to sun direction, in the xy-plane)
rotation_axis = normalize(cross(direction_from_sun_to_earth, [0, 0, 1]))
# Create rotation matrix for 23 degrees south
cos_θ = cos(obliquity_angle)
sin_θ = sin(obliquity_angle)
ux, uy, uz = rotation_axis
# Rodrigues' rotation formula in matrix form
rotation_matrix = [
    cos_θ + ux^2 * (1 - cos_θ)         ux*uy*(1 - cos_θ) - uz*sin_θ    ux*uz*(1 - cos_θ) + uy*sin_θ
    uy*ux*(1 - cos_θ) + uz*sin_θ      cos_θ + uy^2 * (1 - cos_θ)       uy*uz*(1 - cos_θ) - ux*sin_θ
    uz*ux*(1 - cos_θ) - uy*sin_θ      uz*uy*(1 - cos_θ) + ux*sin_θ     cos_θ + uz^2 * (1 - cos_θ)
]
# Apply rotation to get the adjusted sun direction
direction_from_sun_to_earth = rotation_matrix * direction_from_sun_to_earth

sun_light = Makie.DirectionalLight(Colors.WP_A |> RGBf, direction_from_sun_to_earth)
fig = with_theme(theme_dark()) do
    Figure(; figure_padding = 0);
end
ax = GlobeAxis(fig[1, 1]; show_axis = false)
push!(ax.scene.lights, sun_light)
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

satellite_marker = Observable{Point3f}(Point3f(0, 0, 0))
satellite_trajectory = Observable{Vector{Point3f}}(Point3f[(0,0,0), (1, 1, 1)])
# satellite_mesh = load("/Users/anshul/Downloads/ImageToStl.com_CubeSat+-+1+RU+Generic/CubeSat - 1 RU Generic.obj")

satellite_marker_plt = scatter!(
    ax.scene,
    @lift([$satellite_marker]);
    # marker = satellite_mesh,
)

satellite_trajectory_spec = lift(satellite_trajectory) do trajectory
    [Makie.SpecApi.Lines(
        trajectory;
        color = RGBAf.((Makie.wong_colors()[2],), LinRange(0, 1, length(trajectory)))
    )]
end
satellite_trajectory_plt = plotlist!(ax.scene, satellite_trajectory_spec)

record(fig, "earth_rotations.mp4", LinRange(0.01, 5, 2000); framerate = 60, update = false) do t_days
    t_seconds = t_days * 86400
    # Note: we assume the propagator is in the TEME reference frame but the earth is of course in ecef
    # and the sun is in another reference frame.
    # To get all of these together requires some juggling.
    # At some point, GeoMakie's globe axis should just handle all of this
    # and we specify some target datum that it goes to.
    
    # Acquire current GMST so we can adjust our ECEF frame of Earth to the TEME frame (inertial, so Earth rotates within the frame)
    current_gmst = jd_to_gmst(t_days + SpacePowerWorkshop.jd₀)
    # Rotate the earth to reflect the current GMST
    rotate!(earth_plt, Vec3f(0, 0, 1), -current_gmst)
    # Update the satellite marker and trajectory
    current_satellite_position_idx = round(Int, t_seconds) + 1
    satellite_marker[] = SpacePowerWorkshop.sat_pos[current_satellite_position_idx]
    satellite_trajectory[] = SpacePowerWorkshop.sat_pos[max(1, current_satellite_position_idx - 180 * 30):current_satellite_position_idx]
    

end


tmax = 1.0
indmax = 86400 # 1 day in seconds
