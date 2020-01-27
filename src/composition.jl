""" We store the info about the cloud composition in these types.  
The general procedure to fetch the cloud composition at agiven location has two
steps: first we `probe` the cloud at a given location at find the local 
properties relevant to the cloud coposition model (for example this may be the 
droplet radius).  Then we call functions `density`, `radius`, `mie_*` to 
obtain the particular parameter required for some computation.  The reason to
split this procedure into several stages is to avoid repetitive computations
to obtain e.g. the droplet radius.

A cloud-composition type must implement these methods:

 - probe: Return an object that stores the local composition state.  For example, it may return a Float64 with the
    radius of droplets at a given point.  It may also return nothing if we do not need to use any local info.
 - density: droplet density.
 - mie_qext: extinction parameter.
 - mie_g: asymmetry parameter.
 - mie_ω0: Single-scattering albedo.
"""

abstract type CloudComposition end

# Linear regression
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
linreg1(x, y) = x \ y

# Fast computation of x^(-3/4)
power34(x) = 1. / sqrt(x) / sqrt(sqrt(x))


""" A type for a fixed-composition cloud.  It may be initiated with the standard constructor if you know
g, ω0, Qext or with fixedradius if you want to solve the Mie problem for a fixed radius.
"""
struct Fixed <: CloudComposition
    n::Float64
    radius::Float64
    g::Float64
    ω0::Float64
    Qext::Float64
end

probe(comp::Fixed, r::Point) = nothing
density(comp::Fixed, r::Point, ::Nothing) = comp.n
radius(comp::Fixed, r::Point, ::Nothing) = comp.radius
mie_qext(comp::Fixed, r::Point, ::Nothing) = comp.Qext
mie_g(comp::Fixed, r::Point, ::Nothing) = comp.g
mie_ω0(comp::Fixed, r::Point, ::Nothing) = comp.ω0

"""
  Create a FixedComp composition type for a fixed radius and droplet density.
"""
function fixednr(λ, n, radius)
    g, ω0, Qext = MieParams.compute(radius, λ)
    Fixed(n, radius, g, ω0, Qext)
end


struct VariableNR{T,M} <: CloudComposition
    nrfetch::T
    miefit::M
end

probe(comp::VariableNR, r::Point) = radius(comp.nrfetch, r)
density(comp::VariableNR, r::Point, _) = density(comp.nrfetch, r)
mie_qext(comp::VariableNR, r::Point, radius) = mie_qext(comp.miefit, radius)
mie_g(comp::VariableNR, r::Point, radius) = mie_g(comp.miefit, radius)
mie_ω0(comp::VariableNR, r::Point, radius) = mie_ω0(comp.miefit, radius)
