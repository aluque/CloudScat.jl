""" Code to compute Mie scattering parameters based on MieScatter julia
package 

https://github.com/dronir/MieScatter.jl

"""

module MieParams
using FastGaussQuadrature

# Unfortunately the julia package manager does not really work for dependencies
# on unregistered packages.  And it seems that nothing really works for
# sub-packages so I could also not use them.  This was the only half-satisfactory
# way.
include("MieScatter/src/MieScatter.jl")
using .MieScatter

include("constants.jl")
const co = constants

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

end
