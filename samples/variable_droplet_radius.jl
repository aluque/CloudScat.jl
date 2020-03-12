""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("variable_droplet_radius.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia variable_droplet_radius.jl


- - -
2020 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat

const co = CloudScat.constants

# We specify a variable composition through a new type that allows CloudScat
# to dispatch on it whenever it has to access e.g. the radius or density.
struct Composition end

# For the droplet size we will set up a parabolic profile with a maximum at zmax
# where the radius is rmax.
const zmax, rmax = 9 * co.kilo, 20 * co.micro

# The minimum radius is reached at zmin and has a value of rmin.
const zmin, rmin = 12 * co.kilo, 10 * co.micro

# In this case we consider a constant density of scatterers.  One may use
# a more complex method in CloudScat.density below.
const nscat = 100 * co.centi^-3

# Now we overload the CloudScat.radius and CloudScat.density functions to
# dispacth on the type that we have just defined.
function CloudScat.radius(::Composition, r::Point)
    rmax - (rmax - rmin) * (r[3] - zmax)^2 / (zmin - zmax)^2
end
CloudScat.density(::Composition, r::Point) = nscat
    
    
function run()
    params = Params(
        # NUmber of simulated photons
        N = 400000,

        # The wavelength of the scattered light
        λ = 337 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [0, 0, 9 * co.kilo],
        source_b = [0, 0, 9 * co.kilo])
    
    # World geometry
    # Describe the cloud geometry here.  It may consist in combinations of
    # geometrical figures (cylinders, spheres, cones, semi-planes).
    # See cloud_gemoetry.jl for a more detailed description
    #
    # A cylinder is defined as Cylinder(zbottom, ztop, xcenter, ycenter, radius)
    c1 = Cylinder(6 * co.kilo,  12 * co.kilo, 0, 0, 200 * co.kilo)
    cloud = c1
    
    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Cylinder(7 * co.kilo, 60 * co.kilo, 0, 0, 200 * co.kilo)
    
    # We also have to provide a max. collision rate for the null collision
    # method.  The only requirement is that the actual rate never exceeds this
    # but the code is more efficient if the bound is tight.  One option
    # is to provide the highest values of droplet radius and density but
    # if these peak at different places we incur some inefficiency.
    # It works for this case however.
    composition = VariableNR(params.λ, Composition(), rmax, nscat)
                             
    # Define the full simulation world
    world = World(cloud, domain, composition)

    # Observers record a time curve and an image of received light.  You can
    # set up as many observers as you wish but they have a significant impact
    # on performance.
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

    basename = splitext(splitdir(@__FILE__)[2])[1]

    # The saveto option specifies the output file.  Leave it as it is to base
    # the filename on the name of this file specifying the parameters.
    CloudScat.main(params, world, observers,
                   saveto="$(basename).h5")
end

# Check if the code has been 'included' or run from the shell.  In the latter
# case, run the simulation.
isinteractive() || run()
