# Script to invert the size distribution data from Bodega Bay Marine Laboratory

using Dates
using DifferentialMobilityAnalyzers
using Underscores
using Lazy
using RegularizationTools
using Gadfly
using DataFrames
using CSV
using Statistics
using LinearAlgebra

""" df_to_SizeDistribution(df)
    Takes the CSV size distribution format materialized as a DataFrame and converts
    to a vector of SizeDistributions as defined in DifferentialMobilityAnalyzers.jl
"""
function df_to_SizeDistribution(df)
    Dlow = df[1,5:end] |> Vector
    Dp   = df[2,5:end] |> Vector
    Dup  = df[3,5:end] |> Vector
    De   = [Dup[1:end];Dlow[end]]
    ΔlnD = log.(Dup./Dlow)
    map(4:length(df[:,1])) do i
       N = df[i, 5:end] |> Vector
       SizeDistribution([],De, Dp, ΔlnD, N./ΔlnD, N, :response)
    end
end

""" setup_DMA()
    Create setup of DMA as defined in DifferentialMobilityAnalyzers.jl
    Returns Λ, δ as well as cleaned up matrix and number of bins
"""
function setup_DMA()
    lpm = 1.66e-5
    qsa,qsh = 1.3lpm, 5lpm                       
    t,p = 292.15, 1e5                            
    r₁,r₂,l = 9.37e-3,1.961e-2,0.44369           
    leff, m, bins = 73.0, 6, 120 
    DMAtype, polarity = :cylindrical, :-                                

    Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,leff,polarity,m,DMAtype)
    v0 = @>> dtoz(Λ, 12e-9)  ztov(Λ)
    δ = setupSMPS(Λ, v0, 10000, bins, 1)
    A = δ.𝐀[:,:]
    
    Λ, δ, A, bins
end

""" export_to_csv(timestamp, N)
    Create NC State CSV file format for storing size distribution data
    Returns a DataFrame for saving.
"""
function export_to_csv(timestamp, RH, N)
    ts = timestamp
    De = 𝕣[1].De
    Δi = 1
    n = length(ts[1:Δi:end-Δi])

    for i = 1:length(ts)-1
        if ((ts[i+1] - ts[i]) > Dates.Minute(20)) 
            N[i] = zeros(120) 
        end
    end

    yy = mapreduce(vcat, 1:length(ts)) do i
        mmdf = DataFrame(
            Entry = "SMPS Scan", 
            timeISO8601 = ts[i], 
            timeInt64 = Dates.value(ts[i]),
            RH = RH[i]
        )
        map(j -> mmdf[!, Symbol("DMA$j")] = [N[i][j]], 1:length(𝕣[1].Dp))
        mmdf
    end

    @_ filter((sum(_[5:end]) < 10000), eachrow(yy)) |> DataFrame
end

# Define base methods; do not use invert function from Regularizationtools for speed.
Lₖx₀B(ψ, b, x₀, lb, ub) = @> ψ solve(b, x₀, lb, ub) getfield(:x)
LₖB(ψ, b, lb, ub) = @> ψ solve(b, lb, ub) getfield(:x)

df = CSV.read("../../data/raw/bbmlsmps/bbmlsmps.csv", DataFrame)
timestamp = df[4:end,:timeISO8601]
RH = df[4:end,:RH]
𝕣 = df_to_SizeDistribution(df)

Λ, δ, A, bins = setup_DMA()

ψ₀ = setupRegularizationProblem(A, 0)
ψ₂ = setupRegularizationProblem(A, 2)

x₀ = @_ map(inv(δ.𝐒) * _.N, 𝕣)

N0 = (@_ map(Lₖx₀B(ψ₀, 𝕣[_].N, x₀[_], zeros(bins), ones(bins).*Inf), 1:length(𝕣)));
N2 = (@_ map(LₖB(ψ₂, 𝕣[_].N, zeros(bins), ones(bins).*Inf), 1:length(𝕣)));

df0 = @>> export_to_csv(timestamp, RH, N0) vcat(df[1:3,:])
df0 |> CSV.write("../../data/processed/bbmlsmpsL0x0B.csv")

df2 = @>> export_to_csv(timestamp, RH, N2) vcat(df[1:3,:])
df2 |> CSV.write("../../data/processed/bbmlsmpsL2B.csv")