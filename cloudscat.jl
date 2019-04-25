""" Light scattering by clouds

    Alejandro Luque, IAA-CSIC, 2019
"""
module CloudScat

using BenchmarkTools
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

@with_kw mutable struct Params @deftype Float64
    # Cloud limits
    cloud_base = 7 * co.kilo
    cloud_top  = 15 * co.kilo 

    # Location of the source
    source_altitude = 8 * co.kilo
    
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
    nground = 7.767e25

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

    # Maximum number of collisions per photon
    max_iter::Int64 = 1000000
    
    # Domain limit: photons above this height are discarded
    domain_top = 50 * co.kilo

    
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

# Description of the photon population
struct Population
    n::Int64                           # Size of the population
    r::Vector{SVector{3, Float64}}     # Positions
    μ::Vector{SVector{3, Float64}}     # Directions
    t::Vector{Float64}                 # Time
    isactive::Vector{Bool}             # Is active?
end

Population(n::Int64) = Population(n,
                                  Vector{SVector{3, Float64}}(undef, n),
                                  Vector{SVector{3, Float64}}(undef, n),
                                  Array{Float64, 1}(undef, n),
                                  Array{Bool, 1}(undef, n))

struct Observer
    r::SVector{3, Float64}   # Location of the observer
    δt::Float64              # Time-resolution of the detector
    δμ::Float64              # Spatial resoution of the camera (in cosine space)
    μmax::Float64            # Field-of-view in cosine space (i.e. max. observable mu)
    obs::Vector{Float64}     # Observations
    img::Array{Float64, 2}   # Image array
end

# For the dynamical dispatch we define types of collisions but they
# currently don't do anything
abstract type ScatteringType end 
struct Mie <: ScatteringType end
struct Rayleigh <: ScatteringType end
struct Null <: ScatteringType end


# Henyey-Greenstein phase function
phg(g, μ) = (1 / 4π) * (1 - g^2) / (1 + g^2 - 2 * g * μ)^(3/2)

# Sampling from the HG phase function .
# Here ξ in a random variable uniformly distributed between 0 and 1
μhg(g, ξ) = (1 + g^2 - ((1 - g^2) / (1 + g * (2 * ξ - 1)))^2) / (2 * g)

# Rayleigh phase function
pr(μ) = (3 / 16π) * (1 + μ^2)

# Hack to allow interchanging both phase functions
pr(g, μ) = pr(μ)


# Sampling from the Rayleigh phase function.  The PDF notes are wrong.
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


function main(args)
    Crayons.force_color(true)

    infile = args[1]
    outfile = splitext(infile)[1] * ".h5"
    
    input = YAML.load_file(infile)

    paramdict = input["parameters"]
    
    params = Params(;Dict(Symbol(k)=>v for (k, v) in paramdict)...)

    println(BOLD(GREEN_FG(format("[{}] cloudscat.jl {}", Dates.now(), infile))))
    println()
    for field in fieldnames(Params)
        printfmtln("{:<45s}: {}",
                   BOLD(YELLOW_FG(String(field))), getfield(params, field))
    end
    println()

    ## Init all observers
    observers = Vector{Observer}()
        
    for obsdata in input["observers"]
        obs = Observer(obsdata)
        push!(observers, obs)
    end

    ## Init all photons
    p = Population(params.N)
    initphotons!(p, params)
    
    run!(p, observers, params)
    save(outfile, p, observers, params)
    
end


"""
Save the population and observations into a h5 file.
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
            
            g["t", args...] = [obs.δt * (i - 0.5) for i in 1:length(obs.obs)]
            g["timeline", args...] = obs.obs
            g["image", args...] = obs.img
            
        end
    end
    println()
    println(BOLD(GREEN_FG(format("[{}] Output written in {}",
                               Dates.now(), fname))))
    println()
end


""" 
Initialize all photons according to the parameters
"""
function initphotons!(p::Population, params::Params)
    @unpack source_altitude = params

    for i in 1:p.n
        p.r[i] = @SVector [0, 0, source_altitude]
        p.μ[i] = randsphere()
        p.t[i] = 0.0
        p.isactive[i] = true
    end
end


"""
Run the MC simulation on a photon population and a collection of observers.
"""
function run!(p::Population, observers::Vector{Observer}, params::Params)    
    @unpack max_iter, N = params
    
    prog = Progress(N, 5)

    for it in 1:max_iter
        actives = iterate!(p, observers, params)
        
        if actives == 0
            ProgressMeter.finish!(prog)
            break
        end
        if it % 100 == 0
            update!(prog, N - actives,
                    showvalues=[(:iterations, it),
                                (:particles, actives)])
        end

        # Old style progress:
        # if it % 1000 == 0            
        #     println(BLUE_FG(format("[{:<23}]", Dates.now())),
        #             "    iterations: $it; particles: $actives") 
        # end
    end

end


"""
Iterate over all particles in the population p and advance them, including
their eventual observation.
"""
function iterate!(p::Population, observers::Vector{Observer}, params::Params)
    # Active particles    
    c::Int64 = 0
    @inbounds for i in 1:p.n
        p.isactive[i] || continue
        c += 1

        # Propagate
        l = travel(p.r[i][3], p.μ[i][3], params)
        
        p.r[i] = p.r[i] + p.μ[i] * l
        p.t[i] += l / co.c
        #@show p.r[i]
        
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
    c
end


""" 
Samples points from the unitary sphere. 
"""
function randsphere()
    ϕ = 2π * rand()
    sinϕ, cosϕ = sincos(ϕ)
    
    u = 2 * rand() - 1
    v = sqrt(1 - u^2)

    @SVector [v * cosϕ, v * sinϕ, u]
end

"""
Propagate a single particle given its position and direction.
"""
@fastmath function travel(z, μz, params::Params)
    @unpack cloud_top, ν_cloud, Λ_cloud, ν_above, Λ_above = params

    # Select the optical depth to travel
    τ = -log(rand())
    
    if incloud(z, params) && μz > 0
        # τ required to exit the cloud
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
Choose the type of scattering event, depending on the altitude z. 
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
    
    # Now find time of obervation if the particle manages to escape
    tobs = t + sobs / co.c

    # Update the observation curve
    ind = 1 + Int64(fld(tobs, o.δt))
    if ind <= length(o.obs)
        o.obs[ind] += f / o.δt
    end

    # Now update the image at the given pixels
    px = 1 + Int64(fld(μobs[1] + o.μmax, o.δμ))
    py = 1 + Int64(fld(μobs[2] + o.μmax, o.δμ))
    
    if 0 < px <= size(o.img)[1] && 0 < py <= size(o.img)[2]
        # Note that we are not dividing here by the solid angle subtended
        # by the pixel.
        o.img[px, py] += f
    end
end

observeone!(o::Observer, ::Type{Mie}, r, μ, t, params::Params) =
    observeone!(o, phg, r, μ, t, params)
observeone!(o::Observer, ::Type{Rayleigh}, r, μ, t, params::Params) =
    observeone!(o, pr, r, μ, t, params)
function observeone!(o::Observer, ::Type{Null}, r, μ, t, params::Params) end



""" 
Given a unitary vector μ, sample a scattering angle from the HG
distribution and finds a new vector that forms that angle with μ,
which is then overwritten.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, ::Type{Mie}, params)
    @unpack ω₀, g = params

    # Return false if the particle is absorbed
    rand() < ω₀ || return μ, false

    ϕ = 2π * rand()
    cosθ = μhg(g, rand())

    turn!(μ, cosθ, ϕ), true
end


""" 
Given a unitary vector μ, sample a scattering angle from the Rayleigh
distribution and finds a new vector that forms that angle with μ,
which is then overwritten.

NB: This function breaks down if μ is directed exactly along z.
"""
@inline @fastmath function scatterone(μ, ::Type{Rayleigh}, params)
    ϕ = 2π * rand()
    cosθ = μr(rand())

    turn!(μ, cosθ, ϕ), true
end

scatterone(μ, ::Type{Null}, params) = (μ, true)


""" 
Deviate the direction of a particle

Change the unitary vector μ to a new one deviated with inclination θ and azimuth
ϕ.  For performance, the passed parameter is cosθ instead of θ.
"""
@inline @fastmath function turn!(μ, cosθ, ϕ)
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
Reorder particles to have all active particle at the initial positions in
the list.  This improves performance but very little so is not invoked now.

"""
function pack!(p::Population)
    tail = p.n
    for i in 1:p.n
        if i >= tail
            break
        end

        if !p.isactive[i]
            while !p.isactive[tail] && tail > i
                tail -= 1
            end
            if tail != i
                # To keep track of deceased photons we should exchange
                
                p.μ[:, i] .= p.μ[:, tail]
                p.r[:, i] .= p.r[:, tail]
                p.isactive[i] = true
                p.isactive[tail] = false
            end
        end
    end
    n
end

function Observer(obsdata::Dict{Any, Any})
    μmax = sin(deg2rad(obsdata["fov"]))
    δμ = 2 * μmax / obsdata["pixels"]
    
    obs = Observer(SVector(obsdata["shift"], 0,
                           obsdata["altitude"]),
                   obsdata["tsample"],
                   δμ,
                   μmax,
                   zeros(obsdata["nsamples"]),
                   zeros(obsdata["pixels"], obsdata["pixels"]))
end

# For static compilation. See PackageCompiler.jl docs
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    main()
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

end
