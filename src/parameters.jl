export Params, init_params
using REPL


"""
    Params

Structure to contain all simulation parameters.
"""
@kwdef struct Params{T, S}
    " Min z for observations"
    zmin::T = 9 * co.kilo
    
    "Initial segment point in default source"
    source_a::SVector{3, T} = @SVector([-2 * co.kilo, 2 * co.kilo,
                                              13 * co.kilo])

    "End segment point for in default source"
    source_b::SVector{3, T} = @SVector([-2 * co.kilo, 2 * co.kilo,
                                              13 * co.kilo])
                                             
    "Photon source"
    source::S = SegmentImpulsiveSource(source_a, source_b)
    
    "Initial number of photons"
    N::Int64 = 100000

    "Wavelength"
    λ::T = 337 * co.nano
    
    "Minimum weight: below, particles are subjected to Russian roulette"
    weight_min::T = 0.01
    
    "Cross-section for Rayleigh scattering"
    σray::T = Ray.σ(λ)
    
    "Maximum number of collisions per photon"
    max_iter::Int64 = Int(1e9)

    "Air density at ground level (z = 0)"
    nair::T = co.nair

    "Rayleigh scattering rate at ground level"
    νray_ground::T = σray * nair
    
    "Atmospheric scale height"
    H::T = 7.2 * co.kilo
    
    "Max. Rayleigh collision rate"
    νraymax::T = νray_ground * exp(-zmin / H)
    
    # @assert λ > 0
    # @assert H > 0
    # @assert nair > 0
end


function Base.show(io::IO, params::Params; cutoff = 40)
    print(io, "{")
    println(io)

    fwidth = maximum([length(String(field)) for field in fieldnames(Params)]) + 1
    vwidth = maximum([min(cutoff, length(repr(getfield(params, field))))
                      for field in fieldnames(Params)]) + 1
    
    for field in fieldnames(Params)
        doc = string(REPL.fielddoc(Base.Docs.Binding(@__MODULE__, :Params), field))
        # Remove the Default: added by Parameters
        iof = IOBuffer()
        show(iof, getfield(params, field))
        s = String(take!(iof))
        if length(s) > cutoff - 3
            s = s[begin:cutoff - 3] * "..."
        end
        println(io, "    # ", strip(doc))
        printfmtln(io, "    {:<$(fwidth)s}= {:<$(vwidth)s}",
                   String(field), s, doc)
        println(io)
    end
    
    println(io)
    print(io, "}")
end        
