""" 
This is a sample input file for the CloudScat code.  To use it run julia
and type 

julia> include("cloud_geometry.jl"); run()

The output file will be based on this file's name so to run different simulations
you can just copy this file, change the configuration and run it as above with
the new file name.  You can also run from the command line as e.g.

> julia cloud_geometry.jl


- - -
2019 Alejandro Luque IAA-CSIC (aluque@iaa.es)

"""

using CloudScat
using StaticArrays
using LinearAlgebra
using CoordinateTransformations

const co = CloudScat.constants

function run()
    params = Params(
        # NUmber of simulated photons
        N = 1e8,

        # The wavelength of the scattered light
        λ = 337 * co.nano,

        # The source is defined as a straight line running from source_a
        # to source_b.  Within this line photons are distributed uniformly with
        # initially random, isotropic directions.  To use a point source set
        # source_b = source_a
        source_a = [0, 0, 8 * co.kilo],
        source_b = [0, 0, 13 * co.kilo])
    
    # World geometry
    # Describe the cloud geometry here.
    # The geometry is specified using several elementary shapes
    # that can be deformed and combined.  The geometry primitives are defined
    # in src/geometry.jl
    #
    # ELEMENTARY SHAPES
    # =================
    # 
    # Currently the supported elementary shapes are:
    # - A sphere:
    #   Sphere(xcenter, ycenter, zcenter, radius)
    # - A vertical, finite cylinder:
    #   Cylinder(zbottom, ztop, xcenter, ycenter, radius)
    # - A vertical, finite cone (possibly a double-cone):
    #   Cone(zbottom, ztop, xvertex, yvertex, zvertex, m).
    #   m is the slope of the cone, i.e. the cone surface satisfies
    #   sqrt((x - xv)^2 + (y - yv)^2) = m (z - zv).
    # - A half-space delimited by a horizontal plane:
    #   Plane(z, up)
    #   If up is true, the half-space is above the plane, else it is below the
    #   plane
    #
    # COMBINATIONS
    # ============
    # 
    # - Union of several shapes:
    #   union(shape1, shape2, ...)
    #   Alternatively, (shape1 + shape2 + ...) would also work.
    # - Intersection:
    #   intersect(shape1, shape2, ...)
    # - Difference
    #   shapediff(shape1, shape2)
    #   shape1 - shape2 would also work.
    #
    # TRANSFORMATIONS
    # ===============
    #
    # Any linear transformation in the package CoordinateTransformations is
    # supported.  See https://github.com/FugroRoames/CoordinateTransformations.jl
    # To use these transformations use ]add CoordinateTransformations in the julia
    # prompt.
    # Rotations from the Rotations package (]add CoordinateTransformations)
    # are also supported when embedded in a LinearMap.
    # For performance reasons, use StaticArrays whenever possible.

    # As an example we will build a geometry composed by a wide cylinder topped
    # by a conical "turret" which is in turn topped by an ellipsoid and has
    # a small semi-spherical "hole" on top.
    
    # First a 20-km wide cylinder:
    cyl = Cylinder(7 * co.kilo,  10 * co.kilo, 0, 0, 20 * co.kilo)

    # Now a cone on the axis
    cone = Cone(9 * co.kilo, 15 * co.kilo, 0, 0, 25 * co.kilo,
                0.75)

    # Finally to define the ellipsiod we start with a sphere centered on the
    # origin with radius 1 km.
    sph = Sphere(0, 0, 0, 1 * co.kilo)

    # We deform it to an ellipsoid with axes (12 km, 10 km, 2 km) and move it to
    # (0, 0, 15 km)
    map = AffineMap(SMatrix{3, 3, Float64}(Diagonal([12, 10, 2])),
                    SVector{3, Float64}([0, 0, 15 * co.kilo]))
    ellipsoid = TransformedShape(map, sph)
    
    # Now the hole that we will substract from the whole, slightly shifted in the
    # y direction.
    hole = Sphere(0, 1 * co.kilo, 17 * co.kilo, 3 * co.kilo)
    
    # Finally we combine everything together
    cloud = shapediff(union(cyl, cone, ellipsoid), hole)
    
    # Geometry of the full domain.  Photons outside the domain are
    # immediately discarded.
    domain = Cylinder(7 * co.kilo, 60 * co.kilo, 0, 0, 200 * co.kilo)
    
    # Set a homogeneous cloud droplet density and radius.
    # See variable_droplet_radius.jl for how to deal with inhomogeneous clouds.
    # Radius of the scattering particles
    radius = 10e-6

    # Density of the scattering centers.  This is only used to estimate
    # the mean free path in the null collision method; use the highest
    # value here and set a inhomogeneous density using the function nfunc
    # below
    nscat = 100 * co.centi^-3
    
    composition = Homogeneous(params.λ, nscat, radius)

    # Define the full simulation world
    world = World(cloud, domain, composition)

    # Observers record a time curve and an image of received light.  You can
    # set up as many observers as you wish but they have a significant impact
    # on performance.

    observers = [Observer(
        # Location of the observer
        position = [-200 * co.kilo, -200 * co.kilo, 400 * co.kilo],
        
        # Sampling interval
        tsample = 1e-5,
        
        # Collected samples
        nsamples = 10000,
        
        # fov of the camera (half of the diagonal fov).
        fov = 40,
        
        # camera pixels
        pixels = 1024)]

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
    basename = splitext(splitdir(@__FILE__)[2])[1]

    # The saveto option specifies the output file.  Leave it as it is to base
    # the filename on the name of this file specifying the parameters.
    CloudScat.main(params, world, observers,
                   saveto="$(basename).h5")
end

# Check if the code has been 'included' or run from the shell.  In the latter
# case, run the simulation.
isinteractive() || run()
