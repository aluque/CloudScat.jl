""" Code to compute Mie scattering parameters based on MieScatter julia
package 

https://github.com/dronir/MieScatter.jl

"""

module MieScat

using FastGaussQuadrature
using MieScatter

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

end

