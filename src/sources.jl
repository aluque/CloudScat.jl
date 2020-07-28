## Handling of arbitrary sources
## We define sources by implementing the method newphoton, which returns
## the r (location), μ (direction), t (time) and w (weight) of a sampled photon
## for the source

export PointImpulsiveSource, SegmentImpulsiveSource, PointLastingSource,
    GaussianSegmentImpulsiveSource


abstract type AbstractSource end

struct PointImpulsiveSource <: AbstractSource
    r0::Point
end

struct SegmentImpulsiveSource <: AbstractSource
    a::Point
    b::Point
end

struct GaussianSegmentImpulsiveSource <: AbstractSource
    a::Point
    b::Point
end

struct PointLastingSource <: AbstractSource
    r0::Point
    t1::Float64
    t2::Float64
end


""" 
    newphoton(params)

Initialize a a new photon.
"""
function newphoton(s::PointImpulsiveSource)
    r = s.r0
    μ = randsphere()
    t = 0.0
    w = 1.0

    r, μ, t, w
end

function newphoton(s::SegmentImpulsiveSource)
    ξ = trand()
    r = (1 - ξ) * s.a + ξ * s.b
    μ = randsphere()
    t = 0.0
    w = 1.0

    r, μ, t, w
end

function newphoton(s::PointLastingSource)
    r = s.r0
    μ = randsphere()
    t = randlin(s.t1, s.t2)
    w = 1.0

    r, μ, t, w
end

function newphoton(s::GaussianSegmentImpulsiveSource)
    c = 0.5 .* (s.a + s.b)
    r = c .+ 0.5 * randn() .* (s.b - s.a)

    t = 0.0
    w = 1.0

    r, μ, t, w
end


"""
    centroid(s)

Gives the location of the centroid of the source.  Used to compute the delay
to any observer. 
"""
centroid(s::PointImpulsiveSource) = s.r0
centroid(s::SegmentImpulsiveSource) = 0.5 * (s.a + s.b)
centroid(s::GaussianSegmentImpulsiveSource) = 0.5 * (s.a + s.b)
centroid(s::PointLastingSource) = s.r0



##
## Random sampling functions
##

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
    randlin(x1, x2)

Sample points from a distribution with a pdf that grows linearly up to `x1` then
decays linearly up to `x1` + `x2`
"""
function randlin(x1, x2)
    ξ = trand()

    if ξ < x1 / (x1 + x2)
        return sqrt(ξ * x1 * (x1 + x2))
    else
        ξ = ξ - x1 / (x1 + x2)
        return x1 + x2 - sqrt(x2^2 - ξ * x2 * (x1 + x2))
    end
end
