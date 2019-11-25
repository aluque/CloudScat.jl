# Henyey-Greenstein phase function.
phg(g, μ) = (1 / 4π) * (1 - g^2) / (1 + g^2 - 2 * g * μ)^(3/2)

# Sampling from the HG phase function.
# Here ξ in a random variable uniformly distributed between 0 and 1
μhg(g, ξ) = (1 + g^2 - ((1 - g^2) / (1 + g * (2 * ξ - 1)))^2) / (2 * g)

# Rayleigh phase function
pr(μ) = (3 / 16π) * (1 + μ^2)

# Hack to allow interchanging both phase functions
pr(g, μ) = pr(μ)

# Isotropic scattering
piso(g, μ) = (1 / 4π)

# Sample from the Rayleigh phase function.  The PDF notes are wrong.
function μr(ξ)
    z = 2 * (2ξ - 1)
    w = z + sqrt(z^2 + 1)
    w^(1//3) - w^(-1//3)
end
