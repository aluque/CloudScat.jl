""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("sample.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia sample.jl


- - -
2019 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat

const co = CloudScat.constants

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
        source_a = [0, 0, 10 * co.kilo],
        source_b = [0, 0, 10 * co.kilo])
    
    # World geometry
    # Describe the cloud geometry here.  It may consist in combinations of
    # geometrical figures (cylinders, spheres, cones, semi-planes).
    # See cloud_gemoetry.jl for a more detailed description
    #
    # A cylinder is defined as Cylinder(zbottom, ztop, xcenter, ycenter, radius)
    c1 = Cylinder(7 * co.kilo,  12 * co.kilo, 0, 0, 200 * co.kilo)
    cloud = c1
    
    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Cylinder(7 * co.kilo, 60 * co.kilo, 0, 0, 200 * co.kilo)
    
    # Set a homogeneous cloud droplet density and radius.
    # See variable_droplet_radius.jl for how to deal with inhomogeneous clouds.
    # Radius of the scattering particles
    radius = 10e-6
    
    # Density of the scattering centers.
    nscat = 100 * co.centi^-3

    # The refindex keyword is used to specify a data file with the refraction
    # index table, specified as
    # wavelengh (micrometers)   n(real part)   n(imaginary part)
    # For ice use WarrenBrandt.dat
    composition = Homogeneous(params.λ, nscat, radius, refindex="Hale.dat")
    
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
