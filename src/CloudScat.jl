module CloudScat
using StaticArrays
using QuadGK
using LinearAlgebra
using Interpolations
using Dates
using HDF5
using YAML
using Parameters
using Formatting
using Logging
using Dates
using ProgressMeter
using Random

export World, Observer

include("parameters.jl")
include("geometry.jl")
include("phasefuncs.jl")
include("constants.jl")
include("rayleigh.jl")

const co = constants


# For the dynamical dispatch we define types of collisions but they
# currently don't do anything
abstract type ScatteringType end 
struct Mie <: ScatteringType end
struct Rayleigh <: ScatteringType end
struct Null <: ScatteringType end
struct Isotropic <: ScatteringType end


struct World{Tc,Td,Fn}
    cloud::Tc
    domain::Td
    nfunc::Fn
end


"""
    Population

Description of the photon population.
"""
mutable struct Population
    n::Int64                           # Size of the population
    r::Vector{SVector{3, Float64}}     # Positions
    μ::Vector{SVector{3, Float64}}     # Directions
    t::Vector{Float64}                 # Time
    isactive::Vector{Bool}             # Is active?
end


"""
    Population(n)

Create an uninitialized population structure with space to hold up to `n`
photons.
"""
Population(n::Int64) = Population(n,
                                  Vector{SVector{3, Float64}}(undef, n),
                                  Vector{SVector{3, Float64}}(undef, n),
                                  Array{Float64, 1}(undef, n),
                                  Array{Bool, 1}(undef, n))


"""
    Observer

Description of an observer.  An observer is located a some position `r` and
collects photons, building 

 * a timeline of arriving photons (`obs`) with a time resolution `δt`
 * an image of the arriving directions of the photons (`img`) with cosine-angle
   resolution `δμ` and field-of-view `μmax`. 
   
"""
struct Observer
    r::SVector{3, Float64}   # Location of the observer
    δt::Float64              # Time-resolution of the detector
    δu::Float64              # Spatial resoution of the camera (in cosine space)
    umax::Float64            # Field-of-view in cosine space (i.e. max. observable mu)
    obs::Array{Float64, 2}   # Observations
    img::Array{Float64, 3}   # Image array
end


"""
    Observer(obsdata)

Create a observer structure from a dictionary of observer properties.
"""
function Observer(;shift, altitude, fov, pixels, tsample, nsamples)
    loc = SVector(shift, 0, altitude)
    umax = tan(deg2rad(fov)) / sqrt(2)
    δu = 2 * umax / pixels
    nt = Threads.nthreads()

    obs = Observer(loc,
                   tsample,
                   δu,
                   umax,
                   zeros(nsamples, nt),
                   zeros(pixels, pixels, nt))
end


function main(params::Params, world::World, observers::Vector{Observer};
              saveto::Union{String,Nothing}=nothing)
    function fmt(level, _module, group, id, file, line)
        return (:blue, format("{:<23}:", Dates.now()), "")
    end
    logger = ConsoleLogger(meta_formatter=fmt)

    with_logger(logger) do
        # Print a list of parameters
        @info "Input parameters:" params
        
        @info "Running with $(Threads.nthreads()) thread(s) (JULIA_NUM_THREADS)"
        
        @info "Initializing the photon population"
        p = Population(params.N)
        initphotons!(p, params)

        # Run the MC simulation
        @info "Running the MC simulation"
        run!(p, world, observers, params)
        
        !isnothing(saveto) && save(saveto, p, observers, params)
        p, observers
    end
end


"""
    save(fname, p, observers, params)

Save the population and observations into a h5 file called `fname`.
"""
function save(fname, p::Population, observers::Vector{Observer},
              params::Params)
    
    h5open(fname, "w") do file
        g = g_create(file, "parameters")
        for field in fieldnames(Params)
            attrs(g)[String(field)] = getfield(params, field)
        end

        args = ("shuffle", (),
                "deflate", 3)

        # g = g_create(file, "population")

        # g["r", args...] = p.r
        # g["mu", args...] = p.μ
        # g["t", args...] = p.t
        # g["isactive", args...] = p.isactive

        for (i, obs) in enumerate(observers)
            g = g_create(file, format("obs{:05d}", i))
            attrs(g)["altitude"] = obs.r[3]
            attrs(g)["shift"] = obs.r[1]

            rsource = 0.5 * (params.source_a + params.source_b)
            delay = norm(obs.r - rsource) / co.c
            attrs(g)["delay"] = delay
            
            g["t", args...] = [obs.δt * (i - 0.5) for i in 1:size(obs.obs, 1)]

            # The dropdims(sum(...)) is to sum over each thread
            g["timeline", args...] = dropdims(sum(obs.obs, dims=2), dims=2)
            g["image", args...] = dropdims(sum(obs.img, dims=3), dims=3)
            
        end
    end
    @info "Output written in $fname"
end


""" 
    initphotons!(p, params)

Initialize all photons in a population `p` according to the parameters `params`.
"""
function initphotons!(p::Population, params::Params)
    @unpack source_a, source_b = params

    for i in 1:p.n
        ξ = trand()
        p.r[i] = (1 - ξ) * source_a + ξ * source_b
        p.μ[i] = randsphere()
        p.t[i] = 0.0
        p.isactive[i] = true
    end
end


"""
    run!(p, world, observers, interps, params)

Run the MC simulation on a photon population and a collection of observers.
"""
function run!(p::Population, world::World, observers::Vector{Observer},
              params::Params)    
    @unpack max_iter, N = params
    @unpack min_fill_ratio, min_actives_for_repack = params

    set_zero_subnormals(true)
    
    observeall!(p, world, Isotropic, observers, params)

    prog = Progress(N, 5)
    
    for it in 1:max_iter
        actives = iterate!(p, world, observers, params)
        fill_ratio = actives / p.n
        if fill_ratio < min_fill_ratio && actives > min_actives_for_repack
            repack!(p)
        end
            
        if actives == 0
            ProgressMeter.finish!(prog)
            break
        end
        if it % 100 == 0
            update!(prog, N - actives,
                    showvalues=[(:iterations, it),
                                (:particles, actives),
                                (:fill_ratio, fill_ratio)])
        end
    end
end


const count = Int64[]

"""
    iterate!(p, observers, params)

Iterate over all particles in the population `p` and advance them, including
their eventual observation by `observers`.
"""
function iterate!(p::Population, world::World, observers::Vector{Observer},
                  params::Params)
    # Active particles
    # c::Int64 = 0
    
    count .= 0
    
    Threads.@threads for tid in 1:Threads.nthreads()
        @inbounds for i in getrange(tid, p.n)
            p.isactive[i] || continue
            count[Threads.threadid()] += 1

            # Propagate
            l = travel_null(p.r[i], p.μ[i], world, params)
        
            p.r[i] = p.r[i] + p.μ[i] * l
            p.t[i] += l / co.c
        
            # Check if inside the domain
            if !inside(world.domain, p.r[i])
                p.isactive[i] = false
                continue
            end

            # Choose scattering type
            scat = choosescat(p.r[i], world, params)
        
            # Accumulate observations
            if p.r[i][3] > params.zmin
                for o in observers
                    observeone!(o, world, scat, p.r[i], p.μ[i], p.t[i], params)
                end
            end
            
            # Scatter and check if absorbed
            p.μ[i], p.isactive[i] = scatterone(p.μ[i], scat, params)
        end
    end
    sum(count)
end


"""
    observeone!(o, p, r, μ, t, params)

Check one particle and one observer and add the particle's contribution to
the timeline and image.

NOTE: Absorption is not considered here: it would simply add a factor ω₀
"""
@inline @fastmath function observeone!(o::Observer, w::World, p::Function,
                                       r, μ, t, params::Params)
    @unpack g, N, Qext, νray_ground, H, radius = params
    # Distance to the observer
    sobs = norm(o.r - r)
    
    # Director to the observer
    μobs = (o.r - r) / sobs

    μscat = μ ⋅ μobs

    # Optical depth to the observer
    ν(r::Point) = (inside(w.cloud, r) ? (Qext * π * radius^2 * w.nfunc(r))
                   : zero(Float64))

    τmie = optpath(w.cloud, ν, r, o.r)

    # Optical depth to the observer from Rayleigh scattering
    τray = νray_ground * H / μobs[3] * (exp(-r[3] / H) - exp(-o.r[3] / H))

    τ = τmie + τray
    
    # The probability dens that the photon reaches the cloud top w/o scattering
    # NB: p(μ) gives the scattering probablity per unit solid angle,
    # imagine that the photon is a super-particle representing N real photons. Now
    # - the number of photons dispersed in a solid angle dΩ is then N p(μ) dΩ,
    # - the surface at the observer subtended by this solid angle is s^2 dΩ.
    # So the number of photons per unit surface at the observer is
    # N p(μ) / s^2.  To this we have to add an attenuation factor exp(-s / Λ).
    f = p(g, μscat) * exp(-τ) / (sobs^2 * N)
    
    # Now find time of observation if the particle manages to escape
    tobs = t + sobs / co.c

    tid = Threads.threadid()
    
    # Update the observation curve
    ind = 1 + Int64(fld(tobs, o.δt))
    if ind <= size(o.obs, 1)
        @inbounds o.obs[ind, tid] += f / o.δt
    end 

    # Now update the image at the given pixels
    px = 1 + Int64(fld(μobs[1] / μobs[3] + o.umax, o.δu))
    py = 1 + Int64(fld(μobs[2] / μobs[3] + o.umax, o.δu))
    
    if 0 < px <= size(o.img)[1] && 0 < py <= size(o.img)[2]
        # Note that we are not dividing here by the solid angle subtended
        # by the pixel.
        @inbounds o.img[px, py, tid] += f
    end
end

observeone!(o::Observer, w, ::Type{Mie}, r, μ, t, params::Params) =
    observeone!(o, w, phg, r, μ, t, params)
observeone!(o::Observer, w, ::Type{Rayleigh}, r, μ, t, params::Params) =
    observeone!(o, w, pr, r, μ, t, params)
observeone!(o::Observer, w, ::Type{Isotropic}, r, μ, t, params::Params) =
    observeone!(o, w, piso, r, μ, t, params)
function observeone!(o::Observer, w, ::Type{Null}, r, μ, t, params::Params)
end


function observeall!(p::Population, w::World, scat::Type{<:ScatteringType},
                     observers::Vector{Observer}, params::Params)
    Threads.@threads for tid in 1:Threads.nthreads()
        @inbounds for i in getrange(tid, p.n)
            p.isactive[i] || continue
            if p.r[i][3] > params.zmin
                for o in observers
                    observeone!(o, w, scat, p.r[i], p.μ[i], p.t[i], params)
                end
            end
        end
    end
end


""" Compute a range appropriate for this thread. 
    Obtained from https://stackoverflow.com/questions/52593588/julia-why-doesnt-shared-memory-multi-threading-give-me-a-speedup
"""
function getrange(tid, n)
    nt = Threads.nthreads()
    d, r = divrem(n, nt)
    from = (tid - 1) * d + min(r, tid - 1) + 1
    to = from + d - 1 + (tid <= r ? 1 : 0)
    from:to
end


""" 
    randsphere()

Sample points from the unitary sphere. 
"""
function randsphere()
    ϕ = 2π * trand()
    sinϕ, cosϕ = sincos(ϕ)
    
    u = 2 * trand() - 1
    v = sqrt(1 - u^2)

    @SVector [v * cosϕ, v * sinϕ, u]
end


"""
    travel(r, μ, params)

Sample the travel distance a single particle given its position `r` and 
its direction `μ`.  This function considers the interfaces in the cloud geometry
and limits the particle movement up to the next interface.  It also considers
different collision rates in and out of the cloud.

DEPRECATED: DOES NOT WORK
"""
@fastmath function travel(r, μ, world::World, params::Params; ϵ=1e-3)
    @unpack νmax, νraymax = params

    # Select the optical depth to travel
    τ = -log(trand())

    # To check whether we use collision rates in or outside the cloud
    # we test a point slightly ahead of the particle in its path.  This is
    # to allow cases where the particle is exactly at an interface.
    rp = r + (ϵ / νmax) * μ

    # Still to do: use local collision rates.
    ν = inside(world.cloud, rp) ? νmax : νraymax

    l = τ / ν
    rend = r + l * μ

    minimizer = MinAccumulator(1.0)
    interfaces!(minimizer, world.cloud, r, rend)
    l * minimizer.value
end


"""
    travel_null(r, μ, params)

Sample the travel distance a single particle given its position `r` and 
its direction `μ`.
"""
@fastmath function travel_null(r, μ, world::World, params::Params; ϵ=1e-3)
    @unpack νmax = params

    # Select the optical depth to travel
    τ = -log(trand())
    τ / νmax
end



""" 
    choosescat(z, params)

Choose the type of scattering event, depending on the altitude `z`. 
"""
@inline function choosescat(r, world::World, params::Params)
    @unpack νmax, Qext, nair, H, σray, νray_ground, radius = params
    inside(world.domain, r) || return Null

    νMie = (inside(world.cloud, r) ?
            Qext * π * radius^2 * world.nfunc(r) :
            Float64(0))
    νRay = νray_ground * exp(-r[3] / H)
    
    ξ = trand() * νmax

    if ξ <= νMie
        return Mie
    elseif ξ <= (νMie + νRay)
        return Rayleigh
    else
        return Null
    end
end


""" 
Given a unitary vector μ, sample a scattering angle from the HG
distribution and finds a new vector that forms that angle with μ,
which is returned.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, ::Type{Mie}, params)
    @unpack ω₀, g = params

    # Return false if the particle is absorbed
    trand() < ω₀ || return μ, false

    ϕ = 2π * trand()
    cosθ = μhg(g, trand())

    turn(μ, cosθ, ϕ), true
end


""" 
Given a unitary vector μ, sample a scattering angle from the Rayleigh
distribution and finds a new vector that forms that angle with μ,
which is returned.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, ::Type{Rayleigh}, params)
    ϕ = 2π * trand()
    cosθ = μr(trand())

    turn(μ, cosθ, ϕ), true
end


scatterone(μ, ::Type{Null}, params) = (μ, true)


""" 
    turn(μ, cosθ, ϕ)

Find a unitary vector deviated with inclination θ and azimuth
ϕ with respect to μ.  For performance, the passed parameter is cosθ instead of θ.
"""
@inline @fastmath function turn(μ, cosθ, ϕ)
    sinϕ, cosϕ = sincos(ϕ)
    
    sinθ = sqrt(1 - cosθ^2)

    s = sqrt(1 - μ[3]^2)

    b(μ1, μ2, μ3) = μ1 * μ3 * cosϕ + μ2 * sinϕ 

    μ1x = sinθ * b(μ[1], -μ[2], μ[3]) / s + μ[1] * cosθ
    μ1y = sinθ * b(μ[2],  μ[1], μ[3]) / s + μ[2] * cosθ
    μ1z = -s * sinθ * cosϕ + μ[3] * cosθ

    @SVector [μ1x, μ1y, μ1z]
end


"""
    pack!(p)

Reorders the particles in the population `p` to have all active particle 
at the initial positions in the list.

"""
function repack!(p::Population)
    start = time()
    
    # new positions
    k = zeros(Int64, p.n)
    c = 0
    for i in 1:p.n
        if p.isactive[i]
            c += 1
            k[c] = i
        end
    end

    for i in 1:c
        p.μ[i] = p.μ[k[i]]
        p.r[i] = p.r[k[i]]
        p.t[i] = p.t[k[i]] 
        p.isactive[i] = true
    end

    @debug "\u1b[0KParticles repackaged (took $(1000 * (time() - start)) ms)"
    p.n = c
end


const INTERFACE_LISTS = Vector{Float64}[]

@noinline function init_interface_list()
    # Allocate the list on the thread's own heap lazily
    # instead of the master thread heap to minimize memory conflict.
    tid = Threads.threadid()
    INTERFACE_LISTS[tid] = Float64[0.0]
    INTERFACE_LISTS[tid]
end

@inline function get_interface_list()
    @inbounds begin
        tid = Threads.threadid()
        isassigned(INTERFACE_LISTS, tid) || return init_interface_list()
        resize!(INTERFACE_LISTS[tid], 1)
        INTERFACE_LISTS[tid]
    end
end

""" 
Compute the optical path between points `a` and `b` in a given geometry
`geom` and with a position-dependent collision rate `ν`.
"""
function optpath(geom, ν::Function, a::Point, b::Point)
    v = ListAccumulator(get_interface_list())

    interfaces!(v, geom, a, b)
    sort!(v.list)

    push!(v.list, 1.0)
    
    f(s) = ν((1 - s) * a + s * b)
    
    L = norm(a - b)
    # path = L * quadgk(f, 0., v..., 1.0,
    #                   atol=1e-3, rtol=1e-2, order=3)[1]
    path = L * quadgauss(f, v.list)
    path
end


# Pre-computed order-3 quadrature
const quadrule = (@SVector([-0.7745966692414834,
                            0.0,
                            0.7745966692414834]),
                  @SVector([0.5555555555555556,
                            0.8888888888888888,
                            0.5555555555555556]))
@fastmath function quadgauss(f::Function, x)
    xq, wq = quadrule
    x1 = x[1]

    I = zero(x1)
    for ix in 2:length(x)
        h = 0.5 * (x[ix] - x1)
        Ip = zero(x1)
        
        for i in 1:length(wq)
            @inbounds Ip += wq[i] * f(x1 + (1 + xq[i]) * h)
        end
        x1 = x[ix]
        I += h * Ip            
    end
    I
end


""" 
Build a function that returns a fixed number `ν0` inside the geometry 
`geom` and zero otherwise.
"""
function constinside(geom, ν0)
    function ν(r::Point)
        inside(geom, r) ? ν0 : zero(ν0)
    end
end


const RNG = MersenneTwister[]

@noinline function init_thread_rng()
    # Allocate the random number generator on the thread's own heap lazily
    # instead of the master thread heap to minimize memory conflict.
    RNG[Threads.threadid()] = MersenneTwister(Threads.threadid())
end


@inline function trand()
    @inbounds begin
        tid = Threads.threadid()
        isassigned(RNG, tid) || init_thread_rng()
        return rand(RNG[tid])
    end
end


function __init__()
    nth = Threads.nthreads()
    resize!(RNG, nth)
    resize!(INTERFACE_LISTS, nth)
    resize!(count, nth)
    
    # for tid in 1:nth
    #     RNG[tid] = MersenneTwister(tid)
    # end
end


end # module
