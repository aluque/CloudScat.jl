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

abstract type AbstractComposition end

# Linear regression
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
linreg1(x, y) = x \ y

# Fast computation of x^(-3/4)
power34(x) = 1. / sqrt(x) / sqrt(sqrt(x))


""" A type for a fixed-composition cloud.  It may be initiated with the standard constructor if you know
g, ω0, Qext or with fixedradius if you want to solve the Mie problem for a fixed radius.
"""
struct Homogeneous <: AbstractComposition
    n::Float64
    radius::Float64
    g::Float64
    ω0::Float64
    Qext::Float64
    νmiemax::Float64

    function Homogeneous(n, radius, g, ω0, Qext)
        new(n, radius, g, ω0, Qext, Qext * π * radius^2 * n)
    end
end

probe(comp::Homogeneous, r::Point) = nothing
density(comp::Homogeneous, r::Point, ::Nothing) = comp.n
radius(comp::Homogeneous, r::Point, ::Nothing) = comp.radius
mie_qext(comp::Homogeneous, r::Point, ::Nothing) = comp.Qext
mie_g(comp::Homogeneous, r::Point, ::Nothing) = comp.g
mie_ω0(comp::Homogeneous, r::Point, ::Nothing) = comp.ω0
νmiemax(comp::Homogeneous) = comp.νmiemax

"""
  Create a Homogeneous composition instance for a given wavelength `λ`and
  fixed droplet density `n` and droplet radius `r`.
"""
function Homogeneous(λ, n, radius; refindex="Hale.dat")
    interp = m_interpolator(refindex)
    m = interp(λ)
    
    g, ω0, Qext = MieParams.compute(radius, λ, m)
    @info "Mie parameters" g ω0 Qext
    
    Homogeneous(n, radius, g, ω0, Qext)
end


struct VariableNR{T,M} <: AbstractComposition
    nrfetch::T

    # The assymetry parameter is fit as g = mr / (r + r0)
    g0::Float64
    r0::Float64

    # The single-scattering albedo is fit as 1 - ω0 = ar
    a::Float64

    # The Extinction coefficient is fit as Qext = 2 + cr^-3/4
    c::Float64

    νmiemax::Float64
end


"""
    m = VariableNR(λ)

Compute the fit parameters to compute g, ω0 and Qext for an arbitrary radius
in the range 1 μm < r < 100 μm.
"""
function VariableNR(λ, nrfetch, νmiemax; refindex="Hale.dat")
    g0, r0, a, c = miefitparams(λ, refindex)
    VariableNR(nrfetch, g0, r0, a, c, νmiemax)
end

function VariableNR(λ, nrfetch, rmax, nmax; datarefindex="Hale.dat")
    g0, r0, a, c = miefitparams(λ, refindex)
    Qext_max = 2. + comp.c * power34(rmax)
    νmiemax = Qext_max * π * rmax^2 * nmax
end

function miefitparams(λ, refindex)
    interp = m_interpolator(refindex)
    m = interp(λ)

    r = 1 * co.micro:1 * co.micro:100 * co.micro    
    g, ω0, Qext = compute(r, λ, m; order=10000)

    # Linear regression always more robust so we fit r/g ~ r/m + r0 / m
    β = linreg(r, r ./ g)
    g0, r0 = 1 / β[2], β[1] / β[2]

    a = linreg1(r, 1 .- ω0)[1]
    c = linreg1(power34.(r), Qext .- 2)[1]

    (g0, r0, a, c)
end    

probe(comp::VariableNR, r::Point) = radius(comp.nrfetch, r)
density(comp::VariableNR, r::Point, _) = density(comp.nrfetch, r)
radius(comp::VariableNR, r::Point, rad) = rad
mie_qext(comp::VariableNR, r::Point, rad) = 2. + comp.c * power34(rad)
mie_g(comp::VariableNR, r::Point, rad) = m.g0 * rad / (rad + m.r0)
mie_ω0(comp::VariableNR, r::Point, rad) = 1. - m.a * rad
νmiemax(comp::VariableNR) = comp.νmiemax


# Compute a linear interpolator for the refraction index from a datafile
function m_interpolator(fname)
    # If the provided filename does not exist we look for it in our
    # data directory.
    isfile(fname) || (fname = joinpath(@__DIR__, "..", "data", fname))
    
    λ1, n1, k1 = eachcol(readdlm(fname, comments=true))

    # Lambda is in micrometers in the input file.
    λ1 *= co.micro

    interp_n1 = LinearInterpolation(λ1, n1)
    interp_k1 = LinearInterpolation(λ1, k1)
    
    function interp(λ)
        # Here we ensure type stability of the return value
        a::Float64 = interp_n1(λ)
        b::Float64 = interp_k1(λ)
        complex(a, b)
    end
end
