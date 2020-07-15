module CloudScat
using StaticArrays
using LinearAlgebra
using Interpolations
using Dates
using HDF5
using Parameters
using Formatting
using Logging
using Dates
using ProgressMeter
using Random
using Interpolations
using FastGaussQuadrature
using DelimitedFiles

export Params, Point, World, Observer, Homogeneous, VariableNR
export NoAbsorption, StratifiedAbsorption
export time, photo, image, stime

const Point = SVector{3, Float64}

include("constants.jl")
include("parameters.jl")
include("mieparams.jl")
include("sources.jl")
include("composition.jl")
include("absorption.jl")
include("geometry.jl")
#include("phasefuncs.jl")
include("rayleigh.jl")

using .MieParams

const co = constants
const SAMPLES_PATH = normpath(joinpath(@__DIR__, "..", "samples"))
const UTIL_PATH = normpath(joinpath(@__DIR__, "..", "util"))

# For the dynamical dispatch we define types of collisions
abstract type ScatteringType end 
struct Mie <: ScatteringType
    g::Float64
    ω0::Float64
end
struct Rayleigh <: ScatteringType end
struct Null <: ScatteringType end
struct Isotropic <: ScatteringType end

const NonNull = Union{Mie, Rayleigh, Isotropic}

# Phase functions
p(s::Mie, μ) = (1 / 4π) * (1 - s.g^2) / (1 + s.g^2 - 2 * s.g * μ)^(3/2)
p(s::Rayleigh, μ) = (3 / 16π) * (1 + μ^2)
p(s::Isotropic, μ) = (1 / 4π)

# Sampling from the phase functions.  Here ξ is a random number uniformly
# distributed in [0, 1).
sample(s::Mie, ξ) = (1 + s.g^2 - ((1 - s.g^2) / (1 + s.g * (2 * ξ - 1)))^2) / (2 * s.g)
function sample(s::Rayleigh, ξ)
    z = 2 * (2ξ - 1)
    w = z + sqrt(z^2 + 1)
    w^(1//3) - w^(-1//3)
end
    
# Test absorption of the particle
function reweight(s::Mie, w, wmin)
    w1 = s.ω0 * w
    if w1 < wmin
        w1 = trand() < w1 ? 1.0 : 0.0
    end
    w1
end
reweight(s::Rayleigh, w, wmin) = w
reweight(s::Isotropic, w, wmin) = w
reweight(s::Null, w, wmin) = w

# Whether a particle can scatter
scatters(s::Null) = false
scatters(s::Union{Rayleigh, Mie, Isotropic}) = true


struct World{Td <: AbstractShape, Ts <: AbstractSource,
             Tc <: AbstractShape, Tcomp <: AbstractComposition,
             Tabs <: AbstractAbsorber, N}
    domain::Td
    source::Ts
    cloud::Tc
    comp::Tcomp
    absorber::Tabs
    quadrule::Tuple{SVector{N,Float64}, SVector{N,Float64}}

    function World(domain::Td, source::Ts, cloud::Tc, comp::Tcomp, absorb::Tabs,
                   n::Int) where {Td, Ts, Tc, Tcomp, Tabs}
        x, w = gausslegendre(n)
        quadrule = (SVector{n}(x), SVector{n}(w))
        new{Td, Ts, Tc, Tcomp, Tabs, n}(domain, source, cloud, comp, absorb, quadrule)
    end

    World(domain, source, cloud, comp) =
        World(domain, source, cloud, comp, NoAbsorption(), 3)
    World(domain, source, cloud, comp, absorb) =
        World(domain, source, cloud, comp, absorb, 3)
end



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
function Observer(;position, fov, pixels, tsample, nsamples)
    sposition = SVector{3, Float64}(position)
    umax = tan(deg2rad(fov)) / sqrt(2)
    δu = 2 * umax / pixels
    nt = Threads.nthreads()

    obs = Observer(sposition,
                   tsample,
                   δu,
                   umax,
                   zeros(nsamples, nt),
                   zeros(pixels, pixels, nt))
end


times(obs::Observer) = [obs.δt * (i - 0.5) for i in 1:size(obs.obs, 1)]
timeline(obs::Observer) = dropdims(sum(obs.obs, dims=2), dims=2)
image(obs::Observer) = dropdims(sum(obs.img, dims=3), dims=3)

"""
    projectpx(obs, px, H)

  Project pixel coordinates into real space.  

  Finds a pair `[x, y]` of positions that, at an altitude `H` would be projected
  into pixels `px = [xp, yp]` of the camera.  This is useful to simulate 
  observations with the same geometry as an existing observation.
"""
function projectpx(obs::Observer, px::AbstractVector, H)
    @. -(obs.r[3] - H) * (obs.umax - px * obs.δu) + obs.r[[1, 2]]
end


@with_kw struct Results
    time::Vector{Float64}
    photo::Vector{Float64}
    image::Array{Float64, 2}
    delay::Float64
end

time(r::Results) = r.time
stime(r::Results) = time(r) .- r.delay
photo(r::Results) = r.photo
image(r::Results) = r.image


function main(params::Params, world::World, observers::Vector{Observer};
              saveto::Union{String,Nothing}=nothing,
              verbose=true)

    function fmt(level, _module, group, id, file, line)
        return (:blue, format("{:<23}:", Dates.now()), "")
    end
    logger = ConsoleLogger(meta_formatter=fmt)

    with_logger(logger) do
        if verbose
            @info "CloudScat (c) Alejandro Luque IAA-CSIC, 2020"
            @info "Mie solver by Olli Wilkman (https://github.com/dronir/MieScatter.jl)"
            # Print a list of parameters
            @info "Input parameters:" params
            
            @info "Running with $(Threads.nthreads()) thread(s) (JULIA_NUM_THREADS)"
        
            # Run the MC simulation
            @info "Running the MC simulation"
        end
        run!(world, observers, params; verbose=verbose)

        # Build the results object to return; one for each observer
        results = map(observers) do obs
            rsource = 0.5 * (params.source_a + params.source_b)
            Results(time=[obs.δt * (i - 0.5) for i in 1:size(obs.obs, 1)],
                    # The dropdims(sum(...)) is to sum over each thread
                    photo=dropdims(sum(obs.obs, dims=2), dims=2),
                    image=dropdims(sum(obs.img, dims=3), dims=3),
                    delay=norm(obs.r - rsource) / co.c)
        end
        
        !isnothing(saveto) && save(saveto, observers, params)

        results
    end
    observers
end


"""
    save(fname, p, observers, params)

Save the population and observations into a h5 file called `fname`.
"""
function save(fname, world::World, observers::Vector{Observer}, params::Params)
    
    h5open(fname, "w") do file
        g = g_create(file, "parameters")
        for field in fieldnames(Params)
            attrs(g)[String(field)] = getfield(params, field)
        end

        args = ("shuffle", (),
                "deflate", 3)

        for (i, obs) in enumerate(observers)
            g = g_create(file, format("obs{:05d}", i))
            attrs(g)["altitude"] = obs.r[3]
            attrs(g)["shift"] = obs.r[1]

            delay = norm(obs.r - centroid(world.source)) / co.c
            attrs(g)["delay"] = delay
            
            g["t", args...] = times(obs)

            # The dropdims(sum(...)) is to sum over each thread
            g["timeline", args...] = timeline(obs)
            g["image", args...] = image(obs)
            
        end
    end
    @info "Output written in $fname"
end



"""
    run!(p, world, observers, interps, params)

Run the MC simulation on a photon population and a collection of observers.
"""
function run!(world::World, observers::Vector{Observer}, params::Params;
              verbose=true)
    @unpack N = params
    
    # This is to allow quick debug runs without changing an input file.
    N ÷= parse(Int, get(ENV, "CLOUDSCAT_DEBUG_REDUCTION", "1"))
    
    set_zero_subnormals(true)

    if verbose
        prog = Progress(N, 5)
        prog_lock = Threads.SpinLock()
    end
    
    count = Threads.Atomic{Int}(0)
    
    Threads.@threads for tid in 1:Threads.nthreads()
        @inbounds for i in getrange(tid, N)
            r, μ, t, w = newphoton(world.source)
            
            for o in observers
                observeone!(o, world, Isotropic(), r, μ, t, w, params)
            end

            iterate(r, μ, t, w, world, observers, params)            

            if verbose
                lock(prog_lock)
                ProgressMeter.next!(prog)
                unlock(prog_lock)
            end
        end
    end
end


"""
    iterate(r, μ, t, w, world, observers, params)

Iterate over steps for a single particle with starting position `r`,
direction `μ`, time `t` and weight `w`.  The particles lives in `world` and
may be detected by a list of `observers`.  `params` are the simulation parameters.
"""
function iterate(r, μ, t, w, world::World, observers::Vector{Observer},
                  params::Params)
    @unpack max_iter = params
    
    for it in 1:max_iter    
        (w > 0) || break
        
        # Propagate
        l = travel_null(r, μ, world, params)
        
        r1 = r + μ * l
        t += l / co.c
        
        # Check if inside the domain
        if !inside(world.domain, r1)
            w = 0.0
            break
        end

        # Absorption during the propagation
        w *= transmittance(world.absorber, r, r1, μ)
        r = r1

        # Choose scattering type
        scat = choosescat(r, world, params)
        
        # Accumulate observations
        if r[3] > params.zmin
            for o in observers
                observeone!(o, world, scat, r, μ, t, w, params)
            end
        end
        
        # Scatter and check if absorbed
        μ, w = scatterone(μ, w, scat, params)
    end
end


"""
    observeone!(o, w, scat, r, μ, t, weight, params)

Check one particle and one observer and add the particle's contribution to
the timeline and image.

NOTE: Absorption is not considered here: it would simply add a factor ω₀
"""
@inline @fastmath function observeone!(o::Observer, w::World,
                                       s::NonNull,
                                       r, μ, t, weight, params::Params)
    @unpack N, νray_ground, H = params
    # Distance to the observer
    sobs = norm(o.r - r)
    
    # Director to the observer
    μobs = (o.r - r) / sobs

    μscat = μ ⋅ μobs

    # Optical depth to the observer
    ν(r::Point) = miecollrate(r, w, params)

    τmie = optpath(w.cloud, ν, r, o.r, w.quadrule)

    # Optical depth to the observer from Rayleigh scattering
    τray = νray_ground * H / μobs[3] * (exp(-r[3] / H) - exp(-o.r[3] / H))

    τ = τmie + τray
    
    # Transmittance due to background absorption
    trans = transmittance(w.absorber, r, o.r, μobs)
    
    # The probability dens that the photon reaches the cloud top w/o scattering
    # NB: p(μ) gives the scattering probablity per unit solid angle,
    # imagine that the photon is a super-particle representing N real photons. Now
    # - the number of photons dispersed in a solid angle dΩ is then N p(μ) dΩ,
    # - the surface at the observer subtended by this solid angle is s^2 dΩ.
    # So the number of photons per unit surface at the observer is
    # N p(μ) / s^2.  To this we have to add an attenuation factor exp(-s / Λ).
    f = p(s, μscat) * exp(-τ) * trans / (sobs^2 * N)
    
    # Now find time of observation if the particle manages to escape
    tobs = t + sobs / co.c

    tid = Threads.threadid()
    
    # Update the observation curve
    ind = 1 + Int64(fld(tobs, o.δt))
    if ind <= size(o.obs, 1)
        @inbounds o.obs[ind, tid] += weight * f / o.δt
    end 

    # Now update the image at the given pixels
    # The minus are here because to be more user-friendly we invert the
    # image formed in the sensor.
    px = 1 + Int64(fld(-μobs[1] / μobs[3] + o.umax, o.δu))
    py = 1 + Int64(fld(-μobs[2] / μobs[3] + o.umax, o.δu))
    
    if 0 < px <= size(o.img)[1] && 0 < py <= size(o.img)[2]
        # The solid angle subtended by the pixel is roughly
        # δu^2 * (cos θ)^3, where θ is the angle with the vertical.
        # Actually cos θ = μobs[3].  The units would be here
        # photons / sr / source photon.
        @inbounds o.img[px, py, tid] += weight * f / (o.δu^2 * μobs[3]^3)
    end
end

observeone!(o::Observer, w, ::Null, r, μ, t, weight, params::Params) = nothing


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
    travel(r, μ, params)

Sample the travel distance a single particle given its position `r` and 
its direction `μ`.  This function considers the interfaces in the cloud geometry
and limits the particle movement up to the next interface.  It also considers
different collision rates in and out of the cloud.

DEPRECATED: DOES NOT WORK
"""
@fastmath function travel(r, μ, world::World, params::Params; ϵ=1e-3)
    @unpack νraymax = params
    νmax = νraymax + νmiemax(world.comp)
    
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
    @unpack νraymax = params
    νmax = νraymax + νmiemax(world.comp)
    
    # Select the optical depth to travel
    τ = -log(trand())
    τ / νmax
end



""" 
    choosescat(z, params)

Choose the type of scattering event, depending on the altitude `z`. 
"""
@inline function choosescat(r, world::World, params::Params)::ScatteringType
    @unpack νraymax, nair, H, σray, νray_ground = params
    νmax = νraymax + νmiemax(world.comp)

    inside(world.domain, r) || return Null()

    if inside(world.cloud, r)
        loc = probe(world.comp, r)
        Qext = mie_qext(world.comp, r, loc)
        νMie =  Qext * π * radius(world.comp, r, loc)^2 * density(world.comp, r, loc)
    else
        νMie = 0.0
    end
    
    νRay = νray_ground * exp(-r[3] / H)
    
    ξ = trand() * νmax

    if ξ <= νMie
        g = mie_g(world.comp, r, loc)
        ω0 = mie_ω0(world.comp, r, loc)
        return Mie(g, ω0)
    elseif ξ <= (νMie + νRay)
        return Rayleigh()
    else
        return Null()
    end
end


""" 
Given a unitary vector μ, sample a scattering angle from the HG
distribution and finds a new vector that forms that angle with μ,
which is returned.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, w, scat::ScatteringType, params)
    @unpack weight_min = params
    scatters(scat) || (return μ, w)
    
    # Return false if the particle is absorbed
    w = reweight(scat, w, weight_min)

    ϕ = 2π * trand()
    cosθ = sample(scat, trand())

    turn(μ, cosθ, ϕ), w
end


"""
    miecollrate(r, w, params)

Computes Mie collision rate at a given point.
"""
function miecollrate(r::Point, w::World, params::Params)
    inside(w.cloud, r) || return zero(Float64)
    loc = probe(w.comp, r)

    # Is the computation of ω0, g optimized out?
    Qext = mie_qext(w.comp, r, loc)
    Qext * π * radius(w.comp, r, loc)^2 * density(w.comp, r, loc)
end



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
function optpath(geom, ν::Function, a::Point, b::Point, quadrule)
    v = ListAccumulator(get_interface_list())

    interfaces!(v, geom, a, b)
    sort!(v.list)

    push!(v.list, 1.0)
    
    f(s) = ν((1 - s) * a + s * b)
    
    L = norm(a - b)
    # path = L * quadgk(f, 0., v..., 1.0,
    #                   atol=1e-3, rtol=1e-2, order=3)[1]

    path = L * quadgauss(f, v.list, quadrule)
    path
end


@fastmath function quadgauss(f::Function, x, quadrule)
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
    
    # for tid in 1:nth
    #     RNG[tid] = MersenneTwister(tid)
    # end
end


end # module
