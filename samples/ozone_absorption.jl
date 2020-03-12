""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("ozone_absorption.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia ozone_absorption.jl


- - -
2020 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat
using DataFrames
using DelimitedFiles

const co = CloudScat.constants

# This example depends on the ozone profiles used by MODTRAN, which cannot
# be redistributed.  If you have working MODTRAN installation, modify this
# to point to the folder with the profiles.
const MODTRAN_PROFILES = joinpath(@__DIR__, "profiles")

function run()
    profiles = ["MidLat_Summer.dat",
                "MidLat_Winter.dat",
                "SubArc_Summer.dat",
                "SubArc_Winter.dat",
                "Tropical.dat",
                "US_Std_1976.dat"]

    for prof in profiles
        run_profile(prof)
    end
end


function run_profile(fname)
    df = load_profile(joinpath(MODTRAN_PROFILES, fname))
    @info "Profile $(fname) loaded"
    
    # This is the absorption cross-section for 200 nm, where ozone absorption
    # is lowest.  Obtained from
    # L.T Molina and M.J. Molina, "Absolute absorption cross sections of
    # ozone in the 185- to 350-nm wavelength range,"
    # J. Geophys. Res. 91, 14501-14508 (1986).
    σ = 3.145E-19 * co.centi^2

    # The profiles reach only 100 km but the observer is above this altitude
    # We simply add another point with a negligible O3 density.
    z = [df.Z; 1000.0]
    dens = [df.O3; 0.0] 
    T = [df.T; df.T[end]]
   
    absorb = StratifiedAbsorption(σ,
                                  (@. z * co.kilo),
                                  (@. dens * 1e-5 * co.atm / co.k / T))

    params = Params(
        # NUmber of simulated photons
        N = 1000000,

        # The wavelength of the scattered light
        λ = 200 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [0, 0, 10 * co.kilo],
        source_b = [0, 0, 10 * co.kilo],

        # Minimum weight before a particle is discarded
        weight_min = 1e-7)
    
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

    composition = Homogeneous(params.λ, nscat, radius)
    
    # Define the full simulation world
    world = World(cloud, domain, composition, absorb)

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
    profbasename = splitext(fname)[1]

    # The saveto option specifies the output file.  Leave it as it is to base
    # the filename on the name of this file specifying the parameters.
    CloudScat.main(params, world, observers,
                   saveto="$(basename)_$(profbasename).h5")
end

function load_profile(fname)
    local d, colnames
    
    open(fname) do fd
        for i in 1:2
            readline(fd)
        end
        colnames = map(Symbol, (split(readline(fd))))
        readline(fd)

        d = readdlm(fd)
    end
    # From atm * cm / km to m^-3
    # but let's leave this to the code using the dataframe
    # @. d[:, 5:end] = d[:, 5:end] * 1e-5 * co.atm / co.k * T
    
    dict = Dict(col => data for (col, data) in zip(colnames, eachcol(d)))
    df = DataFrame(;dict...)
end


# Check if the code has been 'included' or run from the shell.  In the latter
# case, run the simulation.
isinteractive() || run()
