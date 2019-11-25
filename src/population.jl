"""
    Population

Description of the photon population.
"""
mutable struct Population
    n::Int64                           # Size of the population
    r::Vector{SVector{3, Float64}}     # Positions
    μ::Vector{SVector{3, Float64}}     # Directions
    t::Vector{Float64}                 # Time
    isactive::Vector{Bool}             # Is active?
end


"""
    Population(n)

Create an uninitialized population structure with space to hold up to `n`
photons.
"""
Population(n::Int64) = Population(n,
                                  Vector{SVector{3, Float64}}(undef, n),
                                  Vector{SVector{3, Float64}}(undef, n),
                                  Array{Float64, 1}(undef, n),
                                  Array{Bool, 1}(undef, n))

""" 
    initphotons!(p, params)

Initialize all photons in a population `p` according to the parameters `params`.
"""
function initphotons!(p::Population, params::Params)
    @unpack source_altitude, source_extension = params

    for i in 1:p.n
        z = source_altitude + (rand() - 0.5) * source_extension
        p.r[i] = @SVector [0, 0, z]
        p.μ[i] = randsphere()
        p.t[i] = 0.0
        p.isactive[i] = true
    end
end
