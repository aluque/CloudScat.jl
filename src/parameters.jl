export Params, init_params
using REPL


"""
    Params

Structure to contain all simulation parameters.
"""
@with_kw_noshow struct Params @deftype Float64
    " Min z for observations"
    zmin = 9 * co.kilo
           
    "Initial number of photons"
    N::Int64 = 100000

    "Wavelength"
    λ = 337 * co.nano
    
    "Minimum weight: below, particles are subjected to Russian roulette"
    weight_min = 0.01
    
    "Cross-section for Rayleigh scattering"
    σray = Ray.σ(λ)
    
    "Maximum number of collisions per photon"
    max_iter::Int64 = Int(1e9)

    "Air density at ground level (z = 0)"
    nair = co.nair

    "Rayleigh scattering rate at ground level"
    νray_ground = σray * nair
    
    "Atmospheric scale height"
    H = 7.2 * co.kilo
    
    "Max. Rayleigh collision rate"
    νraymax = νray_ground * exp(-zmin / H)
    
    @assert λ > 0
    @assert H > 0
    @assert nair > 0
end


function Base.show(io::IO, params::Params)
    print(io, "{")
    println(io)

    fwidth = maximum([length(String(field)) for field in fieldnames(Params)]) + 1
    vwidth = maximum([length(repr(getfield(params, field)))
                      for field in fieldnames(Params)]) + 1
    
    for field in fieldnames(Params)
        doc = string(REPL.fielddoc(Base.Docs.Binding(@__MODULE__, :Params), field))
        # Remove the Default: added by Parameters
        res = findlast("Default:", doc)
        doc = doc[1:res[1] - 1]
        
        printfmtln(io, "    {:<$(fwidth)s}= {:<$(vwidth)s} # {}",
                   String(field), repr(getfield(params, field)), doc)
    end
    
    println(io)
    print(io, "}")
end        
