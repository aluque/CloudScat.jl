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

struct Composition end

# Peak altitude and radius
const zmax, rmax = 9 * co.kilo, 20 * co.micro

# Altitude and value of the smaller radius 
const zmin, rmin = 12 * co.kilo, 10 * co.micro

const nscat = 100 * co.centi^-3

# We use a parabolic dependence
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
        source_b = [0, 0, 9 * co.kilo],

        # Radius of the scattering particles
        # When the radius is variable, use the maximum radius.
        # This is used only to compute the maximum collision rate.
        radius = 20e-6,     

        # Density of the scattering centers.  This is only used to estimate
        # the mean free path in the null collision method; use the highest
        # value here and set a inhomogeneous density using the function nfunc
        # below
        nscat = 100 * co.centi^-3)
    
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
    
    miefit = MieFit(params.λ)
    νmiemax = mie_qext(miefit, rmax) * π * rmax^2 * nscat
    
    composition = VariableNR(Composition(), miefit, νmiemax)
                             
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
