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
    p1 = plot(dfl, x = :gf, y = :N, color = :Color, Geom.step,
        Theme(plot_padding=[2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Apparent Growth Factor (-)"),
        Guide.ylabel("Number concentration (cm⁻³)", orientation = :vertical),
        Guide.xticks(ticks = collect(1.2:0.1:2.0)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels1 ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 1.2, xmax = 2.0))

    xlabels2 = collect(1:0.5:3)
    p2 = plot(dfr, x = :gf, y = :Frequency, color = :Color, Geom.step,
        Theme(plot_padding=[2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Growth Factor (-)"),
        Guide.ylabel("Frequency (-)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.yticks(ticks = collect(0:0.2:1)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels2 ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))
	
    hstack(p2,p1)
end

function f01(Dd, Nt, k, seed, gf0)
    Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, k)
    Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
    𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
    gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, k)
    f = @> zeros(k) setindex!(1.0, argmin(abs.(gf .- gf0)))	
    model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 2.5, k));
    𝕣 = model(𝕟ᶜⁿ, f, Dd, gf) 
    dfl1 = DataFrame(gf = 𝕣.Dp ./ (Dd*1e9), N = clean(𝕣.N), Color = ["model" for i in 𝕣.N])
    dfl2 = DataFrame(gf = gf, N = 𝐀*f, Color = ["𝐀₂P<sub>gf</sub>" for i in gf])

    T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
    cr(zˢ, k) = ztod(Λ₁, 1, zˢ) / ztod(Λ₁, k, zˢ)
    DMA₁(𝕟, zˢ, gf) = @_ map(cr(zˢ, _) ⋅ (gfₖ(Λ₁, zˢ, gf, _) ⋅ (T₁(zˢ, _) * 𝕟)), 1:6)
    itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
    DMA₂(𝕟) = δ₂.𝐎 * 𝕟
    
    zˢ = dtoz(Λ₁, Dd);      
    𝕄 = @_ map(itp(_) |> DMA₂, DMA₁(𝕟ᶜⁿ, zˢ, gf[argmin(abs.(gf .- gf0))]))
    dfl3 = DataFrame(gf = 𝕄[1].Dp./(Dd*1e9), N = 𝕄[1].N, Color = "k = 1")
    dfl4 = DataFrame(gf = 𝕄[2].Dp./(Dd*1e9), N = 𝕄[2].N, Color = "k = 2")
    dfl5 = DataFrame(gf = 𝕄[3].Dp./(Dd*1e9), N = 𝕄[3].N, Color = "k = 3")
    dfl = DataFrame(gf = gf, Frequency = f, Color = "P<sub>gf</sub>")
    
    plot_gf_dual([dfl2;dfl3;dfl4;dfl5], dfl) 
end

set_default_plot_size(18cm, 7cm)
@>> f01(100e-9, 3000.0, 60, 1000, 1.6) Gadfly.draw(SVG("figures/f01.svg"))