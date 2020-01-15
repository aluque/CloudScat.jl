""" Code to compute Mie scattering parameters based on MieScatter julia
package 

https://github.com/dronir/MieScatter.jl

"""

module MieParams

using Parameters
using FastGaussQuadrature
using Interpolations
using MieScatter
using DelimitedFiles
using Formatting

include("constants.jl")
const co = constants

# Linear regression
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
linreg1(x, y) = x \ y

# Fast computation of x^(-3/4)
power34(x) = 1. / sqrt(x) / sqrt(sqrt(x))

struct MieFit
    # The assymetry parameter is fit as g = mr / (r + r0)
    g0::Float64
    r0::Float64

    # The single-scattering albedo is fit as 1 - ω0 = ar
    a::Float64

    # The Extinction coefficient is fit as Qext = 2 + cr^-3/4
    c::Float64
end


function Base.show(io::IO, miefit::MieFit)
    print(io, "(")
    for (i, field) in enumerate(fieldnames(MieFit))
        (i > 1) && print(io, ", ")
        printfmt(io, "{} = {:.3e}", String(field), getfield(miefit, field))
    end
    print(io, ")")
end


# A invalid value for default initialization
MieFit() = MieFit(NaN, NaN, NaN, NaN)


"""
    g, ω₀, Qext = compute(r, λ, m; order=1000)

Compute the parameters `g`, `ω₀`, `Qext` for given radius of spherical particle,
`r`, wavelength `λ` and refractive index m (generally a complex number).
"""
function compute(r::Number, λ, m; order=1000)
    μ, w = gausslegendre(order)
    theta = acos.(μ)
    
    S, _, _, _ = compute_mie(size_parameter(r, λ), m, theta)

    # For some reason compute_mie returns Qext and Qback only when we ask
    # for forward and / or backward scattering    
    _, Qsca, Qext, Qback = compute_mie(size_parameter(r, λ), m, [0, π])
    
    ω₀ = Qsca / Qext
    g = sum(μ .* w .* S[:, 1]) / sum(w .* S[:, 1])

    g, ω₀, Qext
end

function compute(r::AbstractVector, λ, m; kw...)
    w = compute.(r, λ, m; kw...)

    g, ω0, Qext = [[itm[i] for itm in w] for i in 1:3]
    g, ω0, Qext
end

# If m is not provided, we take it from interpolation of Hale data.
compute(r, λ; kw...) = compute(r, λ, interp_m(λ); kw...)

# Find the refraction index by interpolation from the Hale data
function interp_m(λ)
    fname = joinpath(@__DIR__, "..", "data", "Hale.dat")
    λ1, n1, k1 = eachcol(readdlm(fname, comments=true))

    # Lambda is in micrometers in the input file.
    λ1 *= co.micro

    interp_n1 = LinearInterpolation(λ1, n1)
    interp_k1 = LinearInterpolation(λ1, k1)

    # Here we ensure type stability of the return value
    a::Float64 = interp_n1(λ)
    b::Float64 = interp_k1(λ)
    complex(a, b)
end

"""
    m = miefit(λ)

Compute the fit parameters to compute g, ω0 and Qext for an arbitrary radius
in the range 1 μm < r < 100 μm.
"""
function miefit(λ)
    r = 1 * co.micro:1 * co.micro:100 * co.micro    
    g, ω0, Qext = compute(r, λ; order=10000)

    # Linear regression always more robust so we fit r/g ~ r/m + r0 / m
    β = linreg(r, r ./ g)
    g0, r0 = 1 / β[2], β[1] / β[2]

    a = linreg1(r, 1 .- ω0)[1]
    c = linreg1(power34.(r), Qext .- 2)[1]
    MieFit(g0, r0, a, c)
end

end
