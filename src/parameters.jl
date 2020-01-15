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
    
    "Radius of the scattering particles (for inhom clouds, an upper bound)"
    radius = 10e-6

    # To allow for a variable particle radius we fit the mie scattering
    # parameters g, ω0 and Qext.  These are the fitting parameters
    "Limit asymmetry parameter for r -> inf"
    g0 = 0.84

    "Fit parameter for g = g0 * r / (r + r0)"
    r0 = 0.0

    "Slope of the dependence 1 - ω0 = a * r"
    a = 0.0

    "Parameter for Qext - 2 = c * r^(-3/4)"
    c = 0.0

    
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
        
    "Max. Mie collision rate."
    νmiemax = 0.0

    "Max. total collision rate"
    νmax = νmiemax + νraymax

    @assert λ > 0
    @assert radius > 0
    @assert nscat > 0
    @assert H > 0
    @assert nair > 0
end


function init_params(;kw...)
    λ = kw[:λ]
    radius = kw[:radius]
    σray = Ray.σ(λ)

    miefit = MieParams.miefit(λ)
    @info "Mie scattering parameters computed for λ=$(λ*1e9) nm" miefit
    @info "Rayleigh cross section computed for λ=$(λ*1e9) nm" σray
    
    Qext_max = 2. + miefit.c * power34(radius)
    
    kwd = Dict(pairs(kw))
    (:g0 in keys(kw)) || (kwd[:g0] = miefit.g0)
    (:r0 in keys(kw)) || (kwd[:r0] = miefit.r0)
    (:a in keys(kw)) || (kwd[:a] = miefit.a)
    (:c in keys(kw)) || (kwd[:c] = miefit.c)
    
    (:σray in keys(kw)) || (kwd[:σray] = σray)
    (:νmiemax in keys(kw)) || (kwd[:νmiemax] = Qext_max * π * radius^2 * kwd[:nscat])
    
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
