using DifferentialMobilityAnalyzers
using Gadfly
using CSV
using Statistics
using StatsBase
using DataFrames
using Underscores
using Colors
using Lazy
using Printf

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

dfr = CSV.read("../../data/raw/bbmlsmps/bbmlsmps.csv", DataFrame)
df0 = CSV.read("../../data/processed/bbmlsmpsL0x0B.csv", DataFrame)
df2 = CSV.read("../../data/processed/bbmlsmpsL2B.csv", DataFrame)

df2a = mapfoldl(vcat, 4:length(df0[!,1])) do i
    ii = df2[!,:timeInt64] .== df0[i,:timeInt64]
    df2[ii,:]
end

df2b = [df2[1:3,:];df2a]

df0a = mapfoldl(vcat, 4:length(df2b[!,1])) do i
    ii = df0[!,:timeInt64] .== df2b[i,:timeInt64]
    df0[ii,:]
end

df0b = [df0[1:3,:];df0a]

dfra = mapfoldl(vcat, 4:length(df2b[!,1])) do i
    ii = dfr[!,:timeInt64] .== df2b[i,:timeInt64]
    dfr[ii,:]
end

dfrb = [df0[1:3,:];dfra]

N0 = @_ map(df0b[_,5:end] |> sum, 4:length(df0b[!,1]))
N2 = @_ map(df2b[_,5:end] |> sum, 4:length(df2b[!,1]))

Dlow = dfr[1,5:end] |> Vector
Dp   = dfr[2,5:end] |> Vector
Dup  = dfr[3,5:end] |> Vector
De   = [Dup[1:end];Dlow[end]]
ΔlnD = log.(Dup./Dlow)

SpectralDensity(N) = N./ΔlnD

j = 8200
mdf0 = DataFrame(Dp = dfrb[2,5:end] |> Vector, N = dfrb[j,5:end] |> Vector, Color = "Response")

Λ, δ, A, bins = setup_DMA()

xx0 = inv(δ.𝐒)*(dfrb[j,5:end] |> Vector)
mdfi = DataFrame(Dp = df0b[2,5:end] |> Vector, N = xx0 |> SpectralDensity, Color = "x₀")
mdf1 = DataFrame(Dp = df0b[2,5:end] |> Vector, N = df0b[j,5:end] |> Vector |> SpectralDensity, Color = "L₀x₀B<sub>[0,∞]</sub>")
mdf2 = DataFrame(Dp = df2b[2,5:end] |> Vector, N = df2b[j,5:end] |> Vector |> SpectralDensity, Color = "L₂B<sub>[0,∞]</sub>")

println(dfrb[j,:timeISO8601])
println(df0b[j,:timeISO8601])
println(df2b[j,:timeISO8601])

function plot_dual(dfl, dfr)
    colors = ["black", "darkred", "gray", "steelblue3"]
    xlabel = log10.([10, 50, 100, 500])
    lfunx = x->ifelse(sum(x .== xlabel) == 1, @sprintf("%i",exp10(x)), "")

    p1 = plot(dfl, x = :Dp, y = :N, Geom.step,
        Theme(default_color = colorant"black", plot_padding=[2mm, 15mm, 2mm, 2mm]), 
        Guide.xlabel("Apparent +1 Mobility Diameter (nm)"),
        Guide.ylabel("Number concentration (cm⁻³)", orientation = :vertical),
        Guide.xticks(ticks = log10.([10:10:100;200;300;400;500;600])),
        Scale.x_log10(labels = lfunx),
        Coord.cartesian(xmin = log10(10), xmax = log10(600)))

    p2 = plot(dfr, x = :Dp, y = :N, color = :Color, Geom.step,
        Theme(plot_padding=[-7mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Mobility Diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm⁻³)", orientation = :vertical),
        Guide.xticks(ticks = log10.([10:10:100;200;300;400;500;600])),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_log10(labels = lfunx),
        Coord.cartesian(xmin = log10(10), xmax = log10(600)))
	
    hstack(p1,p2)
end

p = plot_dual(mdf0, [mdf1;mdf2;mdfi])
set_default_plot_size(18cm, 7cm)
draw(SVG("figures/f05.svg"), p)
