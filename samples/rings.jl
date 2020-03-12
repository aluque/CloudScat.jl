""" 
This is an example of how to extend CloudScat with new geometrical shapes.
Here we will build a cloud geometry consisting in concentric rings.  This
illustrates both the required extensions to define new geometrical shapes
and also how to include shapes defined implicitly by a mathematical function
(this is S = {x: f(x) < 0}).  This involves repeatedly solving the equation 
f(x) = 0 along observation paths and therefore one should expect longer
simulation times.

To use this code run julia and type 

julia> include("rings.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia rings.jl


- - -
2020 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat
using Parameters, Roots
const co = CloudScat.constants

# We need to import these names besides those exported by CloudScat
import CloudScat: interfaces!, inside, accum!, Accumulator

"""
Geometrical shape that encapsulates concentric rings of period `L`, amplitude
`A` and base altitude `z0`.  The top of the cloud is z = z0 + A * cos(2πρ/L),
with ρ = √(x²+y²).
"""
struct ConcentricRings <: Shape
    "Period (m)"
    L::Float64

    "Amplitude (m)"
    A::Float64

    "Base altitude (m)"
    z0::Float64
end


# This function provides a way to find intersections between the shape and
# straight lines between two points a and b.  We must find values of s such that
# (1 - s)a + sb is an interface point (0 < s < 1).  For each of these points
# we must call accum!(v, s).  Here we use the find_zeros function from the julia
# Roots package.  Note that in some cases some roots may not be found so use this
# with care.
function interfaces!(v::Accumulator, rings::ConcentricRings, a::Point, b::Point)
    function f(s)
        @unpack A, z0, L = rings
        p = @. (1 - s) * a + s * b
        return p[3] - (A * cos(2π * sqrt(p[1]^2 + p[2]^2) / L) + z0)
    end

    for s in find_zeros(f, 0.0, 1.0)
        accum!(v, s)
    end
end

# This function just checks whether a point is inside the shape.
function inside(rings::ConcentricRings, p::Point)
    @unpack A, z0, L = rings
    p[3] < (A * cos(2π * sqrt(p[1]^2 + p[2]^2) / L) + z0)
end

#
# Once a new shape is defined (here, ConcentricRings) one can use it and combine
# with other shapes to define a cloud geometry.
#

function run()    
    params = Params(
        # Number of simulated photons
        N = 400000,

        # The wavelength of the scattered light
        λ = 337 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [0, 0, 10 * co.kilo],
        source_b = [0, 0, 10 * co.kilo],

        # Minimum altitude for observations
        zmin = 7 * co.kilo)
    
    # A cylinder is defined as Cylinder(zbottom, ztop, xcenter, ycenter, radius)
    c1 = Cylinder(7 * co.kilo,  15 * co.kilo, 0, 0, 20 * co.kilo)
    
    rings = ConcentricRings(2 * co.kilo, 2 * co.kilo, 10 * co.kilo)

    # We limit the extension of the concentric rings by intersecting the shape
    # with a 20-km-radius cylinder.
    cloud = intersect(c1, rings)
    
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
