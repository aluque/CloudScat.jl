""" 
This is a sample input file for the CloudScat2 code.  To use it run julia
and type 

julia> include("scan_depths.jl")

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia scan_depths.jl


- - -
2019 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat
using StaticArrays
using Printf
using Logging
using Dates
using Formatting

const co = CloudScat.constants

function runone(lmbd, radius_um, h_km)
    params = Params(
        # NUmber of simulated photons
        N = 1e6,

        # The wavelength of the scattered light
        λ = lmbd * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = @SVector([0, 0, h_km * co.kilo]),
        source_b = @SVector([0, 0, h_km * co.kilo]),

        # Radius of the scattering particles
        radius = radius_um * 1e-6,

        zmin = 13 * co.kilo,
        
        # Density of the scattering centers.  This is only used to estimate
        # the mean free path in the null collision method; use the highest
        # value here and set a inhomogeneous density using the function nfunc
        # below
        nscat = 100 * co.centi^-3)
    
    # World geometry
    # Describe the cloud geometry here.  It may consist in combinations of
    # geometrical figures (cylinders, spheres, cones, semi-planes)
    cloud = Cylinder(7 * co.kilo, 15 * co.kilo, 0, 0, 10 * co.kilo)

    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Cylinder(7 * co.kilo, 60 * co.kilo, 0, 0, 100 * co.kilo)
    
    # Set a homogeneous cloud droplet density and radius.
    # See variable_droplet_radius.jl for how to deal with inhomogeneous clouds.
    composition = Fixed(params)
    
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
                    saveto=@sprintf("%s_%d_%02dum_%02dkm.h5",
                                    basename, lmbd, radius_um, h_km))
end

function run()
    for lmbd in [337, 777]
        for radius_um in [10, 20]
            for h_km in 8:14
                runone(lmbd, radius_um, h_km)
            end
        end
    end
end

isinteractive() || run()
