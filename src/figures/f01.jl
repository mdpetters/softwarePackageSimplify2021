using DifferentialMobilityAnalyzers
using Gadfly
using LinearAlgebra
using Printf
using Distributions
using DataFrames
using MLStyle
using CSV
using Lazy
using Colors
using Underscores
using RegularizationTools
    
include("../commonTDMAfunctions.jl")

function plot_gf_dual(dfl, dfr)
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xlabels1 = collect(1.2:0.2:2)
    p1 = plot(dfl, x=:gf, y=:N, color=:Color, Geom.step,
        Theme(plot_padding=[2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Apparent Growth Factor (-)"),
        Guide.ylabel("Concentration (cm⁻³)", orientation=:vertical),
        Guide.xticks(ticks=collect(1.2:0.1:2.0)),
        Guide.colorkey(title=""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels=x -> x in xlabels1 ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin=1.2, xmax=2.0))

    xlabels2 = collect(1:0.5:3)
    p2 = plot(dfr, x=:gf, y=:Frequency, color=:Color, Geom.step,
        Theme(plot_padding=[2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Growth Factor (-)"),
        Guide.ylabel("Probability Density (-)", orientation=:vertical),
        Guide.xticks(ticks=collect(0.8:0.1:2.5)),
        Guide.yticks(),
        Guide.colorkey(title=""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels=x -> x in xlabels2 ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin=0.8, xmax=2.5))
	
    hstack(p2, p1)
end

function f01(Dd, k, gf0)
    Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, k)
	O(k) = mapfoldl(zs -> (δ₂.Ω(Λ₂, δ₂.Z, zs, k) .* δ₂.Tl(Λ₂, δ₂.Z, k))', vcat, δ₂.Z)

    Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
    𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
    gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, k)
	dg = ge[1:end - 1] .- ge[2:end]
    f = @> zeros(k) setindex!(1.0, argmin(abs.(gf .- gf0)))	
    model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 2.5, k));
    𝕣 = model(𝕟ᶜⁿ, f, Dd, gf) 
    dfl1 = DataFrame(gf=𝕣.Dp ./ (Dd * 1e9), N=clean(𝕣.N), Color=["model" for i in 𝕣.N])
    dfl2 = DataFrame(gf=gf, N=𝐀 * f, Color=["𝐁*P<sub>g</sub>" for i in gf])

    T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
	Π(Λ, δ, k) = (@_ map(ztod(Λ, 1, _), dtoz(Λ, k, δ.Dp * 1e-9))) ./ δ.Dp
	DMA₁(𝕟, zˢ, gf) = @_ map(Π(Λ₁, δ₁, _) ⋅ (gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:6)
	itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
	DMA₂(𝕟, k) = O(k) * 𝕟
    
    zˢ = dtoz(Λ₁, Dd)      
	ℕ = DMA₁(𝕟ᶜⁿ, zˢ, gf[argmin(abs.(gf .- gf0))])
	𝕄 = map(k -> (@> itp(ℕ[k]) DMA₂(k)), 1:6)
	dfl3 = DataFrame(gf=𝕄[1].Dp ./ (Dd * 1e9), N=𝕄[1].N, Color="k = 1")
    dfl4 = DataFrame(gf=𝕄[2].Dp ./ (Dd * 1e9), N=𝕄[2].N, Color="k = 2")
    dfl5 = DataFrame(gf=𝕄[3].Dp ./ (Dd * 1e9), N=𝕄[3].N, Color="k = 3")
    dfl = DataFrame(gf=gf, Frequency=f ./ dg, Color="P<sub>g</sub>")
    return plot_gf_dual([dfl2;dfl3;dfl4;dfl5], dfl) 
end

set_default_plot_size(18cm, 7cm)
@>> f01(100e-9, 60, 1.6) Gadfly.draw(SVG("figures/f01.svg"))
