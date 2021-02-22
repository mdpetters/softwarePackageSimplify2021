using DifferentialMobilityAnalyzers
using Gadfly
using LinearAlgebra
using Printf
using Distributions
using DataFrames
using MLStyle
using BenchmarkTools
using CSV
using Lazy
using Colors
using Underscores
using RegularizationTools
	
function plot_gf_dual(dfl, dfr)
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xlabels = collect(1:0.5:3)
    p1 = plot(dfl, x = :gf, y = :N, color = :Color, Geom.step,
        Theme(plot_padding=[2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Apparent Growth Factor (-)"),
        Guide.ylabel("Number concentration (cm⁻³)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))

    p2 = plot(dfr, x = :gf, y = :Frequency, color = :Color, Geom.step,
        Theme(plot_padding=[2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Growth Factor (-)"),
        Guide.ylabel("Frequency (-)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.yticks(ticks = collect(0:0.1:0.5)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))
	
    hstack(p2,p1)
end

include("../commonTDMAfunctions.jl")

Dd = 100e-9
Nt = 3000.0
k = 30
seed = 1000
gf0 = 1.6               

Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, k)
Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, k)
model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 2.5, k));

pop(val,gf0) = @> zeros(k) setindex!(val, argmin(abs.(gf .- gf0)))	
gfs = [1.0, 1.2, 1.6, 2.1]
vals = [0.5,0.15, 0.10, 0.25]

f = @as x mapreduce(i->pop(vals[i], gfs[i]), hcat, 1:4) sum(x; dims = 2) vec(x)
dfr1 = DataFrame(gf = gf, N = 𝐀*f, Color = "Populations")
dfl1 = DataFrame(gf = gf, Frequency = f, Color = "Populations")

Normalize(x) = x./sum(x)
f = @> (0.7*pdf(Normal(1.3,0.07), gf) + pdf(Normal(1.7,0.2), gf)) Normalize
dfr2 = DataFrame(gf = gf, N = 𝐀*f, Color = "Bimodal")
dfl2 = DataFrame(gf = gf, Frequency = f, Color = "Bimodal")

f = @> pdf(truncated(Normal(1.2,0.2) , 1, 17), gf) Normalize
dfr3 = DataFrame(gf = gf, N = 𝐀*f, Color = "Truncated")
dfl3 = DataFrame(gf = gf, Frequency = f, Color = "Truncated")

j = argmin(abs.(gf .- 1.4))
i = argmin(abs.(gf .- 1.8))
f = @> zeros(k) setindex!(ones(j-i+1), i:j) Normalize
dfr4 = DataFrame(gf = gf, N = 𝐀*f, Color = "Uniform")
dfl4 = DataFrame(gf = gf, Frequency = f, Color = "Uniform")

p = plot_gf_dual([dfr1;dfr2;dfr3;dfr4], [dfl1;dfl2;dfl3;dfl4]) 
set_default_plot_size(18cm, 7cm)
Gadfly.draw(SVG("figures/f02.svg"), p)
