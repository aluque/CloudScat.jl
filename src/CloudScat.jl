""" 
cloudscat.jl

Light scattering by clouds

Alejandro Luque, IAA-CSIC, 2019
"""
module CloudScat

using LinearAlgebra
using StaticArrays
using HDF5
using Dates
using Parameters
using YAML
using Crayons
using Crayons.Box
using Formatting
using ProgressMeter

include("constants.jl")
const co = constants

include("Rayleigh.jl")

"""
    Params

Structure to contain all simulation parameters.

Although in the definition many of the parameters take reasonable default
values it's best to set them reading a .yaml file as described in README.md.
"""
@with_kw mutable struct Params @deftype Float64
    # Cloud limits
    cloud_base = 7 * co.kilo
    cloud_top  = 15 * co.kilo 

    # Location of the source
    source_altitude = 8 * co.kilo

    # Vertical extension of the source
    source_extension = 0
    
    # Initial number of photons
    N::Int64 = 10000

    # Extinction efficiency: generally close to 2 for our configuration
    Qext = 2.0
    
    # Assymetry parameter (see e.g. Jursa Handbook of Geophysics and the
    # Space Environment, chap. 18 or Koshak 1994)
    g = 0.84

    # Single scattering albedo (Koshak 1994)
    ω₀ = 0.99996

    # Radius of the scattering particles
    radius = 10e-6

    # Scatering particle density
    nscat = 100 * co.centi^-3

    # Frequency of Mie scattering
    νMie = (Qext * π * radius^2 * nscat)

    # Mean-free-path between Mie scattering events
    ΛMie = 1 / νMie

    #
    # Parameters for Rayleigh scattering
    #
    # Cross section for molecular scattering
    σ = 3e-30

    # Air density at the ground
    nground = 2.6748e25

    # Atmospheric scale height
    H = 7.2e3

    # Frequency of Rayleigh scattering events at the cloud base

    νRay_base = (σ * nground * exp(-cloud_base / H))

    # Frequency of Rayleigh scattering events at the cloud top
    νRay_top = (σ * nground * exp(-cloud_top / H))

    # Max. Frequency of Total scattering inside the cloud
    ν_cloud = νRay_base + νMie

    # Max. Frequency of Total scattering above the cloud
    ν_above = νRay_top
    
    # Mean free paths between collisions (Mie, Rayleigh and Null)
    Λ_cloud = 1 / ν_cloud
    Λ_above = 1 / ν_above

    # Domain limit: photons above this height are discarded
    domain_top = 50 * co.kilo

    # Maximum number of collisions per photon
    max_iter::Int64 = Int(1e9)

    # Minimum fill ratio of the population arrays
    min_fill_ratio = 0.95

    # Maximum number of active particles for repacking
    min_actives_for_repack::Int64 = 1000
    
    @assert Qext > 0
    @assert g > 0
    @assert ω₀ > 0
    @assert radius > 0
    @assert nscat > 0
    @assert 0 <= cloud_base <= cloud_top
    @assert cloud_base <= source_altitude <= cloud_top
    @assert σ >= 0
    @assert nground > 0
    @assert H > 0
    @assert domain_top >= cloud_top
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
    δμ::Float64              # Spatial resoution of the camera (in cosine space)
    μmax::Float64            # Field-of-view in cosine space (i.e. max. observable mu)
    obs::Array{Float64, 2}   # Observations
    img::Array{Float64, 3}   # Image array
end


"""
    Observer(obsdata)

Create a observer structure from a dictionary of observer properties.
"""
function Observer(obsdata::Dict{Any, Any})
    μmax = sin(deg2rad(obsdata["fov"]))
    δμ = 2 * μmax / obsdata["pixels"]
    nt = Threads.nthreads()
    obs = Observer(SVector(obsdata["shift"], 0,
                           obsdata["altitude"]),
                   obsdata["tsample"],
                   δμ,
                   μmax,
                   zeros(obsdata["nsamples"], nt),
                   zeros(obsdata["pixels"], obsdata["pixels"], nt))
end


# For the dynamical dispatch we define types of collisions but they
# currently don't do anything
abstract type ScatteringType end 
struct Mie <: ScatteringType end
struct Rayleigh <: ScatteringType end
struct Isotropic <: ScatteringType end
struct Null <: ScatteringType end


# Henyey-Greenstein phase function.
phg(g, μ) = (1 / 4π) * (1 - g^2) / (1 + g^2 - 2 * g * μ)^(3/2)

# Sampling from the HG phase function.
# Here ξ in a random variable uniformly distributed between 0 and 1
μhg(g, ξ) = (1 + g^2 - ((1 - g^2) / (1 + g * (2 * ξ - 1)))^2) / (2 * g)

# Rayleigh phase function
pr(μ) = (3 / 16π) * (1 + μ^2)

# Hack to allow interchanging both phase functions
pr(g, μ) = pr(μ)

# Hack to use the HG phase function for isotropic scattering
piso(g, μ) = phg(0.0, μ)

# Sample from the Rayleigh phase function.  The PDF notes are wrong.
function μr(ξ)
    z = 2 * (2ξ - 1)
    w = z + sqrt(z^2 + 1)
    w^(1//3) - w^(-1//3)
end


# A function to check if the particle is inside the cloud
incloud(z::Real, params) = (z > params.cloud_base && z < params.cloud_top)
incloud(r::AbstractArray, params) = incloud(r[3], params)

# A function to check if the particle is inside the domain.
indomain(z::Real, params) = (z > params.cloud_base && z < params.domain_top)
indomain(r::AbstractArray, params) = indomain(r[3], params)


function runfromfile(infile)
    # Always use colors
    Crayons.force_color(true)

    # Print greetings
    println(BOLD(GREEN_FG(format("[{}] Cloud-scattering code by A. Luque, IAA-CSIC (aluque@iaa.es)", Dates.now()))))
    println(BOLD(GREEN_FG(format("[{}] cloudscat.jl {}", Dates.now(), infile))))
    println(BOLD(GREEN_FG(format("[{}] Running with {} thread(s) (JULIA_NUM_THREADS)",
                                 Dates.now(), Threads.nthreads()))))

    outfile = splitext(infile)[1] * ".h5"
    params, observers = readinput(infile)

    # Print a list of parameters
    println()
    for field in fieldnames(Params)
        printfmtln("{:<45s}: {}",
                   BOLD(YELLOW_FG(String(field))), getfield(params, field))
    end
    println()

    # Init all photons
    p = Population(params.N)
    initphotons!(p, params)
    
    # Run the MC simulation
    run!(p, observers, params)

    # Save output
    save(outfile, p, observers, params)    
end


"""
    readinput(input)

Read the input parameters and observer configuration from a .yaml file
´infile´.
"""
function readinput(infile)    
    input = YAML.load_file(infile)

    paramdict = input["parameters"]

    # Allow a simple way to include other definitions
    includes = [YAML.load_file(fname) for fname in get(input, "includes", [])]

    for include in includes
        merge!(paramdict, get(include, "parameters", Dict()))
    end
        
    params = Params(;Dict(Symbol(k)=>v for (k, v) in paramdict)...)


    # Init all observers
    observers = Vector{Observer}()
        
    for obsdata in get(input, "observers", [])
        obs = Observer(obsdata)
        push!(observers, obs)
    end

    # Read observers also from files in includes: [...]
    for include in includes
        for obsdata in get(include, "observers", [])
            obs = Observer(obsdata)
            push!(observers, obs)
        end
    end
    
    params, observers
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

            delay = norm(obs.r - [0, 0, params.source_altitude]) / co.c
            attrs(g)["delay"] = delay
            
            g["t", args...] = [obs.δt * (i - 0.5) for i in 1:size(obs.obs, 1)]

            # The dropdims(sum(...)) is to sum over each thread
            g["timeline", args...] = dropdims(sum(obs.obs, dims=2), dims=2)
            g["image", args...] = dropdims(sum(obs.img, dims=3), dims=3)
            
        end
    end
    println()
    println(BOLD(GREEN_FG(format("[{}] Output written in {}",
                               Dates.now(), fname))))
    println()
end


""" 
    initphotons!(p, params)

Initialize all photons in a population `p` according to the parameters `params`.
"""
function initphotons!(p::Population, params::Params)
    @unpack source_altitude, source_extension = params

    for i in 1:p.n
        z = source_altitude + (rand() - 0.5) * source_extension
        p.r[i] = @SVector [0, 0, z]
        p.μ[i] = randsphere()
        p.t[i] = 0.0
        p.isactive[i] = true
    end
end


"""
    run!(p, observers, params)

Run the MC simulation on a photon population and a collection of observers.
"""
function run!(p::Population, observers::Vector{Observer}, params::Params)    
    @unpack max_iter, N = params
    @unpack min_fill_ratio, min_actives_for_repack = params

    observeall!(p, observers, params)
    
    prog = Progress(N, 5)
    
    for it in 1:max_iter
        actives = iterate!(p, observers, params)

        if (actives / p.n) < min_fill_ratio && actives > min_actives_for_repack
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
                                (:fill_ratio, actives / p.n)])
        end
    end
end


"""
    observeall!(p, observers, params)

Calculate the contribution of photons that are not scattered at all.
This function is called only once after the photon population has been 
initialized.
"""
function observeall!(p::Population, observers::Vector{Observer}, params::Params)
    for i in 1:p.n
        # This should not be neccesary but it's best to keep it just in case
        p.isactive[i] || continue

        for o in observers
            observeone!(o, Isotropic, p.r[i], p.μ[i], p.t[i], params)
        end
    end
end


const count = zeros(Int64, Threads.nthreads())
"""
    iterate!(p, observers, params)

Iterate over all particles in the population `p` and advance them, including
their eventual observation by `observers`.
"""
function iterate!(p::Population, observers::Vector{Observer}, params::Params)
    # Active particles
    # c::Int64 = 0
    
    count .= 0
    
    Threads.@threads for tid in 1:Threads.nthreads()
        @inbounds for i in getrange(tid, p.n)
            p.isactive[i] || continue
            count[Threads.threadid()] += 1

            # Propagate
            l = travel(p.r[i][3], p.μ[i][3], params)
        
            p.r[i] = p.r[i] + p.μ[i] * l
            p.t[i] += l / co.c
        
            # Check if inside the domain
            if !indomain(p.r[i], params)
                p.isactive[i] = false
                continue
            end

            # Choose scattering type
            scat = choosescat(p.r[i][3], params)
        
            # Accumulate observations
            for o in observers
                observeone!(o, scat, p.r[i], p.μ[i], p.t[i], params)
            end
            
            # Scatter and check if absorbed
            p.μ[i], p.isactive[i] = scatterone(p.μ[i], scat, params)
        end
    end
    sum(count)
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
    ϕ = 2π * rand()
    sinϕ, cosϕ = sincos(ϕ)
    
    u = 2 * rand() - 1
    v = sqrt(1 - u^2)

    @SVector [v * cosϕ, v * sinϕ, u]
end


"""
    travel(z, μz, params)

Sample the travel distance a single particle given its position `z` and 
z-component of its direction `μz`.
"""
@fastmath function travel(z, μz, params::Params)
    @unpack cloud_top, ν_cloud, Λ_cloud, ν_above, Λ_above = params

    # Select the optical depth to travel
    τ = -log(rand())
    
    if incloud(z, params) && μz > 0
        # Depth required to exit the cloud
        τ0 = ν_cloud * (cloud_top - z) / μz
        if τ < τ0
        # Stay in the cloud
            return Λ_cloud * τ
        else
            return Λ_cloud * τ0 + Λ_above * (τ - τ0)
        end
    end
    
    if incloud(z, params) && μz <= 0
        # Stay in the cloud or leave it from below, in which case the photon
        # is removed.
        return Λ_cloud * τ
    end

    # The photon is outside the cloud.
    if μz > 0
        # Propagate upwards
        return Λ_above * τ
    end

    # The photon is moving downwards and may enter the cloud
    τ0 = ν_above * (cloud_top - z) / μz
    if τ < τ0
        return Λ_above * τ
    else
        return Λ_above * τ0 + Λ_cloud * (τ - τ0)
    end

end


""" 
    choosescat(z, params)

Choose the type of scattering event, depending on the altitude `z`. 
"""
@inline function choosescat(z, params::Params)
    @unpack cloud_base, cloud_top, ν_cloud, ν_above = params
    @unpack νMie, νRay_base, νRay_top, H = params
    indomain(z, params) || return Null
    
    if incloud(z, params)
        ξ = rand() * ν_cloud
        if ξ <= νMie
            return Mie
        elseif ξ <= (νMie + νRay_base * exp(-(z - cloud_base) / H))
            return Rayleigh
        else
            return Null
        end
    else
        ξ = rand() * ν_above
        if ξ < νRay_top * exp(-(z - cloud_top) / H)
            return Rayleigh
        else
            return Null
        end
    end
end


"""
    observeone!(o, p, r, μ, t, params)

Check one particle and one observer and add the particle's contribution to
the timeline and image.

NOTE: Absorption is not considered here: it would simply add a factor ω₀
"""
@inline @fastmath function observeone!(o::Observer, p,
                                       r, μ, t, params::Params)
    @unpack cloud_top, g, N, H, σ, nground, νMie = params
    
    # Distance to the observer
    sobs = norm(o.r - r)
    
    # Director to the observer
    μobs = (o.r - r) / sobs

    μscat = μ ⋅ μobs

    # Optical depth to the observer from Rayleigh scattering
    τ = σ * nground * H / μobs[3] * (exp(-r[3] / H) - exp(-o.r[3] / H))

    if incloud(r, params)
        # Find the distance to the cloud top
        s = sobs * (cloud_top - r[3]) / (o.r[3] - r[3])

        # Depth due to Mie
        τ += s * νMie
    end
    
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
    if ind <= length(o.obs)
        o.obs[ind, tid] += f / o.δt
    end

    # Now update the image at the given pixels
    px = 1 + Int64(fld(μobs[1] + o.μmax, o.δμ))
    py = 1 + Int64(fld(μobs[2] + o.μmax, o.δμ))
    
    if 0 < px <= size(o.img)[1] && 0 < py <= size(o.img)[2]
        # Note that we are not dividing here by the solid angle subtended
        # by the pixel.
        o.img[px, py, tid] += f
    end
end

observeone!(o::Observer, ::Type{Mie}, r, μ, t, params::Params) =
    observeone!(o, phg, r, μ, t, params)
observeone!(o::Observer, ::Type{Rayleigh}, r, μ, t, params::Params) =
    observeone!(o, pr, r, μ, t, params)
observeone!(o::Observer, ::Type{Isotropic}, r, μ, t, params::Params) =
    observeone!(o, pr, r, μ, t, params)
function observeone!(o::Observer, ::Type{Null}, r, μ, t, params::Params) end



""" 
Given a unitary vector μ, sample a scattering angle from the HG
distribution and finds a new vector that forms that angle with μ,
which is returned.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, ::Type{Mie}, params)
    @unpack ω₀, g = params

    # Return false if the particle is absorbed
    rand() < ω₀ || return μ, false

    ϕ = 2π * rand()
    cosθ = μhg(g, rand())

    turn(μ, cosθ, ϕ), true
end


""" 
Given a unitary vector μ, sample a scattering angle from the Rayleigh
distribution and finds a new vector that forms that angle with μ,
which is returned.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, ::Type{Rayleigh}, params)
    ϕ = 2π * rand()
    cosθ = μr(rand())

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


# For static compilation. See PackageCompiler.jl docs
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    runfromfile(ARGS[1])
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    runfromfile(ARGS[1])
end

end
