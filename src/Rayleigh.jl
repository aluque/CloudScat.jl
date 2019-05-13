"""
Compute Rayleigh scattering cross-sections.

Most of the expressions are obtained from 
Bodhaine et al. (1999) J. Atmosph. Ocean. Tech., 16 1854.

If invoked as julia Rayleigh.jl, computes the cross sections for 337 and 777 nm.
From the julia REPL you can compute cross sections for any wavelength as e.g.
```julia-repl
julia> include("Rayleigh.jl")
julia> using Main.Rayleigh
julia> rayleigh(337 * co.nano)
3.4371080894577986e-30
```
Alejandro Luque, IAA-CSIC, 2019
"""
module Rayleigh

export AIR_STP, rayleigh, co

include("constants.jl")
const co = constants

# Depolarization function for different gases.  λ here is μm.
const F = Dict("N2" => λ -> 1.034 + 3.17e-4 / λ^2,
               "O2" => λ -> 1.096 + 1.385e-3 / λ^2 + 1.448e-4 / λ^4,
               "Ar" => λ -> 1.00,
               "CO2" => λ -> 1.15)


# These are densities from ICAO, extracted from Wikipedia.  Note that
# the expressions in Bodhaine assume 300 ppmv CO2 but we consider this
# close enough for our purposes.
const AIR_STP = Dict("N2" => 780840,
                     "O2" => 209476,
                     "Ar" => 9340,
                     "CO2" => 314)

# The density of air in constants.jl is defined at (273.15 K, 101.325 kPa)
# but the functions in Bodhaine et al. 1999 are valid for 288.15 K so we have to
# rescale
nair = co.nair * 273.15 / 288.15


""" 
Compute the refractive index of air.

Takes the wavelength λ in m and returns (n - 1).
Uses the formula from Edlén (1953) referenced from Bodhaine et al. 1999.

"""
function edlen(λ)
    a, b, c = 6432.8, 2949810, 25540
    p, q = 146, 41

    λ_μm = λ / co.micro
    λ2 = 1 / λ_μm^2

    1.0 + 1e-8 * (a + b / (p - λ2) + c / (q - λ2))
end

""" 
Compute the refractive index of air.

Takes the wavelength λ in m and returns (n - 1).
Uses the formula from Peck and Reeder (1953) referenced from Bodhaine et al. 1999.

"""
function peck(λ)
    a, b, c = 8060.51, 2480990, 17455.7
    p, q = 132.274, 39.32957

    λ_μm = λ / co.micro
    λ2 = 1 / λ_μm^2

    1.0 + 1e-8 * (a + b / (p - λ2) + c / (q - λ2))
end


"""
Compute the Rayleigh scattering cross-section.
"""
function rayleigh(λ; composition=AIR_STP)
    n = peck(λ)
    P = sum(ppmv * F[gas](λ / co.micro) for (gas, ppmv) in composition) / 1e6
    
    N = (n^2 - 1)^2 / (n^2 + 2)^2
    24π^3 * N * P / (λ^4 * nair^2) 
end

function main()
    println("# Rayleigh scattering cross-section for 337 nm")
    σ = rayleigh(337 * co.nano)
    println("σ: $σ")

    println()
    println("# Rayleigh scattering cross-section for 777 nm")
    σ = rayleigh(777 * co.nano)
    println("σ: $σ")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# Depolarization functions 
end
