""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("first_hitting.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia first_hitting.jl


- - -
2010 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat
using FastGaussQuadrature

const co = CloudScat.constants

function run()
    params = Params(
        # NUmber of simulated photons
        N = 5000000,

        # The wavelength of the scattered light
        λ = 337 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [0, 0, 10 * co.kilo],
        source_b = [0, 0, 10 * co.kilo],

        # For this example we neglect Rayleigh scattering
        σray = 0.)

    cth = 15 * co.kilo

    # World geometry
    # Describe the cloud geometry here.  It may consist in combinations of
    # geometrical figures (cylinders, spheres, cones, semi-planes).
    # See cloud_gemoetry.jl for a more detailed description
    #
    cloud = Plane(cth, false)
    
    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Plane(60 * co.kilo, false)
        
    # Set a homogeneous cloud droplet density and radius.
    # See variable_droplet_radius.jl for how to deal with inhomogeneous clouds.
    # Radius of the scattering particles
    radius = 10e-6
    
    # Density of the scattering centers.
    nscat = 100 * co.centi^-3

    composition = Homogeneous(params.λ, nscat, radius)
    
    # Define the full simulation world
    world = World(cloud, domain, composition)

    # To capture all exiting photons we set a few observers around the closest
    # point in th cloud top to the source.  If the radiance would be perfectly
    # lambertian this would not be neccesary.
    μ, w = gausslegendre(5)
    
    observers = @. observerat(acos(0.5 * (μ + 1)), cth)

    basename = splitext(splitdir(@__FILE__)[2])[1]

    # The saveto option specifies the output file.  Leave it as it is to base
    # the filename on the name of this file specifying the parameters.
    CloudScat.main(params, world, observers,
                   saveto="$(basename).h5")
end


function observerat(θ, zc)
    R = 400 * co.kilo

    Observer(
        # Location of the observer
        position = [0, R * sin(θ), zc + R * cos(θ)],
        
        # Sampling interval
        tsample = 1e-5,
        
        # Collected samples
        nsamples = 10000,
        
        # fov of the camera (half of the diagonal fov).
        fov = 40,
        
        # camera pixels
        pixels = 1024)
end


# Check if the code has been 'included' or run from the shell.  In the latter
# case, run the simulation.
isinteractive() || run()
