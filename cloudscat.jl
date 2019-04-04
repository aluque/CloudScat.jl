""" Light scattering by clouds

    Alejandro Luque, IAA-CSIC, 2019
"""
module CloudScat

using BenchmarkTools
using LinearAlgebra
using StaticArrays
using HDF5

include("constants.jl")
const co = constants

# Assymetry parameter (see e.g. Jursa Handbook of Geophysics and the
# Space Environment, chap. 18 or Koshak 1994)
const g = 0.84

# Single scattering albedo (Koshak 1994)
const ω₀ = 0.99996

# Radius of the scattering particles
const radius = 10e-6

# Scatering particle density
const nscat = 100 * co.centi^-3

# Mean-free-path between scatterings
# Koshak uses half this but I do not understand why
const Λ = 1 / (π * radius^2 * nscat)

# Maximum number of photons
const N = 1000_000

# Henyey-Greenstein phase function
p(μ) = 0.5 * (1 - g^2) / (1 + g^2 - 2 * g * μ)^(3/2)

# Sampling from the HG phase function (note that this is slightly different than
# in Koshak.  Here ξ in a random variable uniformly distributed between 0 and 1
hgμ(ξ) = (1 + g^2 - ((1 - g^2) / (1 + g * (2 * ξ - 1)))^2) / (2 * g)

# Travel distance from a random variable ξ in (0, 1)
travel(ξ) = -Λ * log(ξ)

# A function to checkk if the particle is inside the cloud.
const cloud_base = 7 * co.kilo
const cloud_top  = 15 * co.kilo 
incloud(r) = (r[3] > cloud_base && r[3] < cloud_top)

# Description of the photon population
struct Population
    r::Array{Float64, 2}     # Positions
    μ::Array{Float64, 2}     # Directions
    t::Vector{Float64}     # Time
    isactive::Vector{Bool} # Is active?
end

struct Observer
    r::SVector{3, Float64}   # Location of the observer
    δt::Float64              # Time-resolution of the detector
    obs::Vector{Float64}     # Observations
end

function main()
    @show nscat
    @show radius
    @show Λ
    @show N
    
    nobs = 10_000
    asim = Observer(SVector(300 * co.kilo, 0, 400 * co.kilo),
                    10 * co.micro,
                    zeros(nobs))

    # Storage space for photon positions and directions and active/inactive
    r = Array{Float64, 2}(undef, 3, N)
    μ = Array{Float64, 2}(undef, 3, N)
    t = Array{Float64, 1}(undef, N)
    isactive = Array{Float64, 1}(undef, N)

    # This is to later allow for reordering
    n = N
    
    for i in 1:n
        r[:, i] .= [0, 0, 8 * co.kilo]
        randsphere!(@view μ[:, i])
        t[i] = 0.0
        isactive[i] = true
    end

    p = Population(r, μ, t, isactive)

    for it in 1:1000000
        actives = iterate!(p, asim, n)
        
        if actives == 0
            break
        end

        if it % 1000 == 0
            @show it, actives
        end
    end

    h5open("/tmp/cloudscat.h5", "w") do file
        write(file, "population/r", p.r)
        write(file, "population/mu", p.μ)
        write(file, "population/t", p.t)
        write(file, "population/isactive", p.isactive)
        write(file, "observer/t",
              [asim.δt * (i - 0.5) for i in 1:length(asim.obs)])
        write(file, "observer/timeline", asim.obs)
    end        
end


"""
Iterate over all particles in the population p and advance them, including
their eventual observation.
"""
function iterate!(p::Population, o::Observer, n::Int64)
    # Active particles
    c::Int64 = 0
    @inbounds for i in 1:n
        p.isactive[i] || continue
        c += 1

        # Propagate
        l = travel(rand())
        p.r[:, i] .+= p.μ[:, i] .* l
        p.t[i] += l / co.c

        # Check if inside the cloud
        incloud(p.r[:, i]) || (p.isactive[i] = false)

        # Accumulate observation
        observeone!(o, p.r[:, i], p.μ[:, i], p.t[i])

        # Scatter and check if absorbed
        p.isactive[i] = scatterone!(@view p.μ[:, i])
    end
    c
end


""" 
Samples points from the unitary sphere. 
"""
function randsphere!(μ)
    ϕ = 2π * rand()
    sinϕ, cosϕ = sincos(ϕ)
    
    u = 2 * rand() - 1
    v = sqrt(1 - u^2)
    
    μ[1] = v * cosϕ
    μ[2] = v * sinϕ
    μ[3] = u
end

"""
Checks one particle and one observer and adds the particle's contribution to
the timeline.

NOTE: Absorption is not considered here: it would simply add a factor ω₀
"""
function observeone!(o::Observer, r, μ, t)
    # Distance to the observer
    sobs = norm(o.r - r)
    
    # Director to the observer
    μobs = (o.r - r) ./ sobs

    μscat = μ ⋅ μobs
    
    # Find the distance to the cloud top
    s = sobs * (cloud_top - r[3]) / (o.r[3] - r[3])

    # The probability dens that the photon reaches the cloud top w/o scattering
    # NB: p(μ) gives the scattering probablity per unit solid angle so one
    # should not include here the sin(θ) factor.
    f = p(μscat) * exp(-s / Λ) / (4π * sobs^2)

    # Now find time of obervation if the particle manages to escape
    tobs = t + sobs / co.c

    # Update the observation curve
    ind = 1 + Int64(floor(tobs / o.δt))
    if ind <= length(o.obs)
        o.obs[ind] += f
    end
    
end


""" 
Given a unitary vector μ, sample a scattering angle from the HG
distribution and finds a new vector that forms that angle with μ,
which is then overwritten.

NB: This function breaks down if μ is directed exactly along z.
"""
function scatterone!(μ)
    # Return false if the particle is absorbed
    rand() < ω₀ || return false

    ϕ = 2π * rand()
    cosθ = hgμ(rand())

    sinϕ, cosϕ = sincos(ϕ)
    
    sinθ = sqrt(1 - cosθ^2)

    s = sqrt(1 - μ[3]^2)

    b(μ1, μ2, μ3) = μ1 * μ3 * cosϕ + μ2 * sinϕ 

    μ1x = sinθ * b(μ[1], -μ[2], μ[3]) / s + μ[1] * cosθ
    μ1y = sinθ * b(μ[2],  μ[1], μ[3]) / s + μ[2] * cosθ
    μ1z = -s * sinθ * cosϕ + μ[3] * cosθ

    μ[1] = μ1x
    μ[2] = μ1y
    μ[3] = μ1z

    true
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
