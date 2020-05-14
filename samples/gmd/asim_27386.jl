""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("asim_27386.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia asim_27386.jl


- - -
2019 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat
using StaticArrays
using LinearAlgebra
using CoordinateTransformations

const co = CloudScat.constants

# We specify a variable composition through a new type that allows CloudScat
# to dispatch on it whenever it has to access e.g. the radius or density.
struct Composition end

# For the droplet size we will set up a linear profile with a maximum at zmax
# where the radius is rmax.
const zmax, rmax = 15 * co.kilo, 20 * co.micro

# The minimum radius is reached at zmin and has a value of rmin.
const zmin, rmin = 5 * co.kilo, 10 * co.micro


# In this case we consider a constant density of scatterers.  One may use
# a more complex method in CloudScat.density below.
const nscat = 100 * co.centi^-3

# Now we overload the CloudScat.radius and CloudScat.density functions to
# dispacth on the type that we have just defined.
function CloudScat.radius(::Composition, r::Point)
    rmin + (r[3] - zmin) / (zmax - zmin) * (rmax - rmin)
end
CloudScat.density(::Composition, r::Point) = nscat


function run()
    observers = [Observer(
        # Location of the observer
        position = [0, 0, 400 * co.kilo],
        
        # Sampling interval
        tsample = 1e-5,
        
        # Collected samples
        nsamples = 10000,
        
        # fov of the camera (half of the diagonal fov).
        fov = 40,
        
        # camera pixels
        pixels = 1024)]


    # We obtain approximate physical coordinates from the pixel coordinates
    # of the observation
    x, y = CloudScat.projectpx(observers[1], [1008, 743], 10 * co.kilo)
    
    params = Params(
        # NUmber of simulated photons
        N = 5e8,

        # The wavelength of the scattered light
        λ = 337 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [x + 1 * co.kilo, y + 2 * co.kilo, 7 * co.kilo],
        source_b = [x + 1 * co.kilo, y + 2 * co.kilo, 7 * co.kilo])
    
    # First a 40-km wide cylinder:
    cyl = Cylinder(5 * co.kilo,  12 * co.kilo, x, y, 15 * co.kilo)

    sph = Sphere(0, 0, 0, 1 * co.kilo)
    
    map0 = AffineMap(SMatrix{3, 3, Float64}(Diagonal([9, 2, 2])),
                     SVector{3, Float64}([x + 2 * co.kilo, y - 4 * co.kilo, 12 * co.kilo]))
    depr = TransformedShape(map0, sph)
    
    # Now the turret
    # turret = Cylinder(9 * co.kilo,  15 * co.kilo, x, y, 5 * co.kilo)
    map1 = AffineMap(SMatrix{3, 3, Float64}(Diagonal([4, 6, 2])),
                     SVector{3, Float64}([x - 1 * co.kilo, y - 1 * co.kilo, 12 * co.kilo]))
    turret1 = TransformedShape(map1, sph)

    map2 = AffineMap(SMatrix{3, 3, Float64}(Diagonal([4, 5, 5])),
                     SVector{3, Float64}([x - 4.0 * co.kilo, y + 3 * co.kilo,
                                          12 * co.kilo]))
    turret2 = TransformedShape(map2, sph)
    
    
    # Finally we combine everything together
    cloud = union(shapediff(cyl, depr), turret1, turret2)
    
    viewpoint = "{" * join(observers[1].r, ", ") * "}"
    viewdir = "{" * join([observers[1].r[1], observers[1].r[2], 0.0],  ", ") * "}"
    println("To visualize the cloud geometry with Mathematica use:")
    println()
    println("R = ", mathematica(cloud))
    println("""RegionPlot3D[R, PlotTheme -> "Classic", Axes -> True, 
        PlotStyle -> Directive[Opacity[0.75], RGBColor[.6, .65, .85], 
            Specularity[0.2]], 
        Lighting -> "Neutral", PlotPoints -> 100, 
        ViewVector -> {$viewpoint, $viewdir}]
    """)
    println()

    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Cylinder(7 * co.kilo, 60 * co.kilo, x, y, 200 * co.kilo)
    
    # We also have to provide a max. collision rate for the null collision
    # method.  The only requirement is that the actual rate never exceeds this
    # but the code is more efficient if the bound is tight.  One option
    # is to provide the highest values of droplet radius and density but
    # if these peak at different places we incur some inefficiency.
    # It works for this case however.
    composition = VariableNR(params.λ, Composition(), rmax, nscat)
    @show composition
    
    # Define the full simulation world
    world = World(cloud, domain, composition)

    # Observers record a time curve and an image of received light.  You can
    # set up as many observers as you wish but they have a significant impact
    # on performance.
    basename = splitext(splitdir(@__FILE__)[2])[1]

    # The saveto option specifies the output file.  Leave it as it is to base
    # the filename on the name of this file specifying the parameters.
    CloudScat.main(params, world, observers,
                   saveto="$(basename).h5")
end

# Check if the code has been 'included' or run from the shell.  In the latter
# case, run the simulation.
isinteractive() || run()
