""" Code to compute Mie scattering parameters based on MieScatter julia
package 

https://github.com/dronir/MieScatter.jl

"""

module MieParams

using FastGaussQuadrature
using Interpolations
using MieScatter
using DelimitedFiles

include("constants.jl")
const co = constants


"""
    g, ω₀, Qext = compute(r, λ, m; order=1000)

Compute the parameters `g`, `ω₀`, `Qext` for given radius of spherical particle,
`r`, wavelength `λ` and refractive index m (generally a complex number).
"""
function compute(r, λ, m; order=1000)
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

# If m is not provided, we take it from interpolation of Hale data.
compute(r, λ; kw...) = compute(r, λ, interp_m(λ); kw...)

function interp_m(λ)
    fname = joinpath(@__DIR__, "..", "data", "Hale.dat")
    λ1, n1, k1 = eachcol(readdlm(fname, comments=true))

    # Lambda is in micrometers in the input file.
    λ1 *= co.micro

    interp_n1 = LinearInterpolation(λ1, n1)
    interp_k1 = LinearInterpolation(λ1, k1)

    complex(interp_n1(λ), interp_k1(λ))
end

end
