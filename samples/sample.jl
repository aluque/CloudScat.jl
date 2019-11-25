""" 
This is a sample input file for the CloudScat2 code.  To use it run julia
and type 

julia> include("sample.jl")

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia sample.jl


- - -
2019 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat2
using StaticArrays

const co = CloudScat2.constants

using Logging
using Dates
using Formatting

function fmt(level, _module, group, id, file, line)
    return (:blue, format("{:<23}:", Dates.now()), "")
end

const logger = ConsoleLogger(meta_formatter=fmt)

function run()
    params = init_params(
        # NUmber of simulated photons
        N = 400000,

        # The wavelength of the scattered light
        λ = 337 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = @SVector([0, 0, 10 * co.kilo]),
        source_b = @SVector([0, 0, 10 * co.kilo]),

        # Radius of the scattering particles
        radius = 10e-6,     

        # Density of the scattering centers.  This is only used to estimate
        # the mean free path in the null collision method; use the highest
        # value here and set a inhomogeneous density using the function nfunc
        # below
        nscat = 100 * co.centi^-3)
    
    # World geometry
    # Describe the cloud geometry here.  It may consist in combinations of
    # geometrical figures (cylinders, spheres, cones, semi-planes)
    c1 = Cylinder(7 * co.kilo,  12 * co.kilo, 0, 0, 200 * co.kilo)
    cloud = c1
    
    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Cylinder(7 * co.kilo, 60 * co.kilo, 0, 0, 200 * co.kilo)
    
    # A function to set inhomogeneous densities of scattering particles
    # Here we just use a homogeneous one.
    nfunc(r::Point) = params.nscat
    
    # Define the full simulation world
    world = World(cloud, domain, nfunc)

    # Observers record a time curve and an image of received light.  You can
    # set up as many observers as you wish but they have a significant impact
    # on performance.

    observers = [Observer(
        # Altitude of the observer
        altitude = 400 * co.kilo,
        
        # Horizontal distance
        shift = 0,
        
        # Sampling interval
        tsample = 1e-5,
        
        # Collected samples
        nsamples = 10000,
        
        # fov of the camera
        fov = 40,
        
        # camera pixels
        pixels = 1024)]

    basename = splitext(splitdir(@__FILE__)[2])[1]

    # The saveto option specifies the output file.  Leave it as it is to base
    # the filename on the name of this file specifying the parameters.
    CloudScat2.main(params, world, observers,
                    saveto="$(basename).h5")
end

with_logger(logger) do
    run()
end