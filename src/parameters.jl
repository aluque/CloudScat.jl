export Params, init_params
using REPL

include("mieparams.jl")


"""
    Params

Structure to contain all simulation parameters.
"""
@with_kw_noshow struct Params @deftype Float64
    " Min z for observations"
    zmin = 9 * co.kilo
    
    # The source is assumed to be a line joining points a and b below
    
    "Location of point a of the source"
    source_a::SVector{3, Float64} = @SVector [-2 * co.kilo, 2 * co.kilo,
                                              13 * co.kilo]

    "Location of point b of the source"
    source_b::SVector{3, Float64} = @SVector [-2 * co.kilo, 2 * co.kilo,
                                              13 * co.kilo]
    
    "Initial number of photons"
    N::Int64 = 100000

    "Wavelength"
    λ = 337 * co.nano
    
    "Extinction efficiency"
    Qext = 2.0
    
    "Asymmetry parameter" #(see e.g. Jursa Handbook of Geophysics and the
    # Space Environment, chap. 18 or Koshak 1994)
    g = 0.84

    "Single scattering albedo"
    ω₀ = 0.99996

    "Radius of the scattering particles"
    radius = 10e-6

    "Scatering particle density"
    nscat = 100 * co.centi^-3

    "Cross-section for Rayleigh scattering"
    σray = 0.0
    
    "Domain limit: photons above this height are discarded"
    domain_top = 50 * co.kilo

    "Maximum number of collisions per photon"
    max_iter::Int64 = Int(1e9)

    "Minimum fill ratio of the population arrays"
    min_fill_ratio = 0.95

    "Maximum number of active particles for repacking"
    min_actives_for_repack::Int64 = 100
    
    "Air density at ground level (z = 0)"
    nair = co.nair

    "Rayleigh scattering rate at ground level"
    νray_ground = σray * nair
    
    "Atmospheric scale height"
    H = 7.2 * co.kilo
    
    "Max. Rayleigh collision rate"
    νraymax = νray_ground * exp(-zmin / H)
        
    "Max. Mie collision rate"
    νmiemax = Qext * π * radius^2 * nscat

    "Max. total collision rate"
    νmax = νmiemax + νraymax

    @assert λ > 0
    @assert Qext > 0
    @assert g > 0
    @assert ω₀ > 0
    @assert radius > 0
    @assert nscat > 0
    @assert H > 0
    @assert nair > 0
end


function init_params(;kw...)
    λ = kw[:λ]
    radius = kw[:radius]
    g, ω₀, Qext = MieParams.compute(radius, λ)
    σray = Ray.σ(λ)
    
    @info "Scattering parameters computed for radius=$(radius*1e6) μm, λ=$(λ*1e9) nm" g ω₀ Qext σray

    kwd = Dict(pairs(kw))
    (:g in keys(kw)) || (kwd[:g] = g)
    (:ω₀ in keys(kw)) || (kwd[:ω₀] = ω₀)
    (:Qext in keys(kw)) || (kwd[:Qext] = Qext)
    (:σray in keys(kw)) || (kwd[:σray] = σray)

    Params(;kwd...)
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
