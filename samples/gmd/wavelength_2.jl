""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("wavelength.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia wavelength.jl


- - -
2010 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat

const co = CloudScat.constants

function run()
    for h in [10 * co.kilo, 12 * co.kilo],
        λ in [337 * co.nano, 777 * co.nano],
        R in [10 * co.micro, 20 * co.micro]
        runone(λ, h, R)
    end
end

function runone(λ, h, R)
    params = Params(
        # Number of simulated photons
        N = 10000000,

        # The wavelength of the scattered light
        λ = λ,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [0, 0, h],
        source_b = [0, 0, h])

    cth = 15 * co.kilo

    # World geometry
    # Describe the cloud geometry here.  It may consist in combinations of
    # geometrical figures (cylinders, spheres, cones, semi-planes).
    # See cloud_gemoetry.jl for a more detailed description
    #
    cloud = Cylinder(5 * co.kilo,  cth, 0, 0, 20 * co.kilo)
    
    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = shapediff(Plane(60 * co.kilo, false), Plane(7 * co.kilo, false))
        
    # Density of the scattering centers.
    nscat = 100 * co.centi^-3

    composition = Homogeneous(params.λ, nscat, R)
    
    # Define the full simulation world
    world = World(cloud, domain, composition)

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
    ofile = ("$(basename)_$(Int(λ/co.nano))nm_" *
             "$(Int(h/co.kilo))km_$(Int(R/co.micro))um.h5")
    CloudScat.main(params, world, observers,
                   saveto=ofile)
end


# Check if the code has been 'included' or run from the shell.  In the latter
# case, run the simulation.
isinteractive() || run()
