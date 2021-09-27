using DifferentialMobilityAnalyzers
using Gadfly
using LinearAlgebra
using Printf
using Distributions
using DataFrames
using MLStyle
using LsqFit
using CSV
using Lazy
using Colors
using Underscores
using Random
using RegularizationTools
	
function plot_gf_dual(dfl, dfr)
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xlabels = collect(1:0.5:3)
    p1 = plot(dfl, x = :gf, y = :N, color = :Color, Geom.step,
        Theme(plot_padding=[0mm, 4mm, 2mm, 2mm]), 
        Guide.xlabel("Apparent Growth Factor (-)"),
        Guide.ylabel("Concentration (cm⁻³)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))

        p2 = plot(dfr, x = :gf, y = :Frequency, color = :Color, Geom.step,
        Theme(plot_padding=[-2mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Growth Factor (-)"),
        Guide.ylabel("Probability Density (-)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))
	
    hstack(p1,p2)
end

include("../commonTDMAfunctions.jl")

Qcpc = 1.0 # Flow in LPM
Dd = 100e-9
Nt = 3000.0
k = 60
seed = 1000
gf0 = 1.6               
gf1 = 1.6               

Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, k)
Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, k)
model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 5.0, k));

dg = ge[1:end-1] .- ge[2:end]
f = @> zeros(k) setindex!(1.0, argmin(abs.(gf .- gf0)))	

N0 = 𝐀*f
N1 = poisson_noise(Qcpc, N0; seed = 707, t = 2.0)
N1 = poisson_noise(Qcpc, N0; seed = 714, t = 2.0)

x₀ = N1./sum(N1)
lb, ub = zeros(k), ones(k)
xλ1 = invert(𝐀, N1, Lₖx₀B(2, x₀, lb, ub))
e1 = @> sqrt.(sum((xλ1 .- f).^2.0)./k) round(digits = 3)
xλ2 = invert(𝐀, N1, LₖDₓB(0, 0.001, lb, ub))
e2 = @> sqrt.(sum((xλ2 .- f).^2.0)./k) round(digits = 3)
 
fun(_, p) = (model(𝕟ᶜⁿ, p[1], Dd, p[2])).N
fit = curve_fit(fun, N1, N1, [1, 1.2], lower = zeros(2))
Ax = fit.param
kk = argmin(abs.(Ax[2] .- gf))
xλ3 = @> zeros(length(gf)) setindex!(Ax[1], kk)
e3 = @> sqrt.(sum((xλ3 .- f).^2.0)./k) round(digits = 3)

dfl1 = DataFrame(gf = gf, Frequency = f ./ dg, Color = "Truth")
dfl2 = DataFrame(gf = gf, Frequency = xλ1 ./ dg, Color = "L<sub>2</sub>x<sub>0</sub>B<sub>[0,1]</sub>, $(e1)")
dfl3 = DataFrame(gf = gf, Frequency = xλ2 ./ dg, Color = "L<sub>0</sub>D<sub>1e-3</sub>B<sub>[0,1]</sub>, $(e2)")
dfl4 = DataFrame(gf = gf, Frequency = xλ3 ./ dg, Color = "LSQ<sub>1</sub>, $(e3)")

dfr1 = DataFrame(gf = gf, N = N1, Color = "Input")
dfr2 = DataFrame(gf = gf, N = 𝐀*xλ1, Color = "𝐁*L<sub>2</sub>x<sub>0</sub>B<sub>[0,1]</sub>")
dfr3 = DataFrame(gf = gf, N = 𝐀*xλ2, Color = "𝐁*L<sub>0</sub>D<sub>1e-3</sub>B<sub>[0,1]")
dfr4 = DataFrame(gf = gf, N = 𝐀*xλ3, Color = "𝐁*LSQ<sub>1</sub>")


p = plot_gf_dual([dfr1;dfr2;dfr3;dfr4], [dfl1;dfl2;dfl3;dfl4]) 
set_default_plot_size(18cm, 7cm)
Gadfly.draw(SVG("figures/f04.svg"), p)
