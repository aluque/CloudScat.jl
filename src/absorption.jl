""" Here we define types to handle molecular absorption. """

abstract type AbstractAbsorber end

struct NoAbsorption <: AbstractAbsorber end

# The only type that a subtype of ContinuumAbsorber must implement is
# transmittance, which gives the transmitance between two points

""" Transmittance between points a and b.  The vector μ is passed for performance,as it is often already computed in the calling environment.  It must be
μ = (b - a) / |b - a|. """
transmittance(::NoAbsorption, a, b, μ) = 1.0


# Now we define a more useful absorption for a species with a single absorption
# cross-section and a stratified density
struct StratifiedAbsorption{I<:AbstractInterpolation} <: AbstractAbsorber
    "The cross-section of the absorption."
    σ::Float64

    "An interpolator for the cumulative density of the absorbing species."
    cuminterp::I

    # Using the cumulative density allows us to avoid integrals.  We define
    # cuminterp(z) = ∫(from z0 to z) n(z') dz' so the integral from z1 to z2 is
    # cuminterp(z2) - cuminterp(z1).
end


""" Initialize a StratifiedAbsorption with the given cross-section and 
    number density evaluated at the values of z.
"""
function StratifiedAbsorption(σ, z::AbstractVector, n::AbstractVector)
    cumn = cumtrapz(z, n)
    cuminterp = LinearInterpolation(z, cumn)
    StratifiedAbsorption(σ, cuminterp)
end


function transmittance(str::StratifiedAbsorption, a, b, μ)
    vertint = str.cuminterp(b[3]) - str.cuminterp(a[3])
    
    if μ[3] != 0
        return @fastmath exp(-str.σ * abs(vertint / μ[3]))
    else
        # This is the unlikely situation where the ray propagates exaclty
        # horizontally.  This complicates things a bit.
        # Since we only store the cumulative integral we obtain the density
        # by evaluating the derivative.
        n = Interpolations.gradient(a[3])[1]
        return @fastmath exp(-str.σ * n * norm(a - b))
    end
end
    

struct MultipleAbsorption{T<:Tuple} <: AbstractAbsorber
    absorbers::T
end


# See geometry.jl for a similar construction and why we use generated functions
@generated function transmittance(u::MultipleAbsorption{T}, a, b, μ) where {T}
    L = fieldcount(T)
    out = quote end

    push!(out.args, :(total = 1.0))
    
    for i in 1:L
        push!(out.args,
              quote
              total *= transmittance(u.absorbers[$i], a, b, μ)
              end
              )
    end
    push!(out.args, :(return total))

    out
end


""" Cumulative trapezoidal integration of a function with values `f` evaluated 
    at `x`.  The result is stored in `r`.
"""
function cumtrapz!(r, x, f)
    @assert length(x) == length(f)
    @assert length(r) == length(r)
    cum = zero(eltype(r))
    
    for i in eachindex(x)
        r[i] = cum
        if i != lastindex(x)
            cum += 0.5 * (f[i + 1] + f[i]) * (x[i + 1] - x[i])
        end
    end
    r
end

""" Like `cumtrapz!` but allocates the output vector. """
function cumtrapz(x, f)
    r = similar(f)
    cumtrapz!(r, x, f)
end

