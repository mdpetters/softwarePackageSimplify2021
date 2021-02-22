# Invert all HTDMA data. This script is slow since it inverts a multi-week dataset

using DifferentialMobilityAnalyzers
using Gadfly
using LinearAlgebra
using Printf
using Distributions
using DataFrames
using Random
using LsqFit
using MLStyle
using Optim
using CSV
using RegularizationTools
using Lazy
using MLStyle
using Colors
using Underscores
using Dates
using DataStructures
using Interpolations
using NumericIO
using ProgressMeter

include("sgpLoaders.jl")
include("../commonTDMAfunctions.jl")

function plot_gf_dual(dfl, dfr)
    set_default_plot_size(18cm, 7cm)
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xlabels = collect(1:0.5:3)
    p1 = plot(dfl, x = :gf, y = :N, color = :Color, Geom.step,
		Guide.xlabel("Growth Factor (-)"),
		Guide.ylabel("Number (cm-3)", orientation = :vertical),
		Guide.xticks(ticks = collect(0.8:0.1:2)),
		Scale.color_discrete_manual(colors...),
		Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
		Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits=1)),
		Theme(plot_padding=[2mm,-6mm,2mm,2mm]),
		Coord.cartesian(xmin = 0.8, xmax = 2))

    colors = ["darkred", "steelblue3", "darkgoldenrod"]
    p2 = plot(dfr, x = :gf, y = :Frequency, color = :Color, Geom.step,
		Guide.xlabel("Growth Factor (-)"),
		Guide.ylabel("Frequency (-)", orientation = :vertical),
  		Guide.xticks(ticks = collect(0.8:0.1:2)),
		Scale.color_discrete_manual(colors...),
		Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
		Theme(plot_padding=[6mm,2mm,2mm,2mm]),
		Coord.cartesian(xmin = 0.8, xmax = 2))
	
    hstack(p1,p2)
end

function plot_psd(dfl, dfr)
	colors = ["black", "darkred"]
    set_default_plot_size(18cm, 7cm)
	
	xt = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]
	xlabels = log10.([10, 20, 50, 100, 200, 500])
    
	p1 = plot(dfl, x = :Dp, y = :S, color = :Dist, Geom.step,
		Guide.xlabel("Particle diameter (nm)"),
		Guide.ylabel("Number (cm-3)", orientation = :vertical),
		Scale.color_discrete_manual(colors...),
        Guide.colorkey(; title = "GF Raw"),
        Guide.xticks(ticks = log10.(xt)),
		Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits=1)),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
		Theme(plot_padding=[2mm,0mm,2mm,2mm]),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)))
	
    p2 = plot(dfr, x = :Dp, y = :S, color = :Dist, Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm-3)"),
        Guide.xticks(ticks = log10.(xt)),
        Guide.colorkey(; title = "PSD"),
		Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits=1)),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
		Theme(plot_padding=[5mm,2mm,2mm,2mm]),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)))
	
	hstack(p2, p1)
end

function prepInversion(i)
	Dd = htdma.Dd[i]*1e-9
	n = length(htdma.Dp[:,i])
	δ₂ = setupDMA(Λ₂, dtoz(Λ₂, 2.5*Dd), dtoz(Λ₂, 0.8*Dd), bins)

    j = argmin(abs.(smps.timestamp .- htdma.timestamp[i]))
    𝕟ᶜⁿ = interpolateSizeDistributionOntoδ((smps.𝕟[j],δ₁))

    gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, bins)
    model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 2.5, bins));
    
 	itp1 = interpolate((htdma.gf[:,i],), htdma.N[:,i], Gridded(Linear()))
 	etp = extrapolate(itp1, 0.0) 
 	R = etp(gf)
	
 	tscn = Dates.Time(smps.timestamp[j])
 	tsht = Dates.Time(htdma.timestamp[i])

 	cn_label = ["$(tscn)" for i = 1:length(𝕟ᶜⁿ.Dp)]
 	ht_label = ["$(tsht)" for i = 1:length(δ₂.Dp)]
 	Ddnm = @sprintf("%i", Dd*1e9)
 	Dd_label = ["Dd = $(Ddnm) nm" for i = 1:2]
 	dfp1 = DataFrame(Dp = 𝕟ᶜⁿ.Dp, S = 𝕟ᶜⁿ.S, Dist = cn_label)
 	dfp2a = DataFrame(Dp = δ₂.Dp, S = R, Dist = ht_label)
 	dfp2b = DataFrame(Dp = [Dd, Dd]*1e9, S = [0, maximum(R)], Dist = Dd_label)
    
    p = plot_psd([dfp2a;dfp2b], dfp1)

    𝐀, gf, ge, R, 𝕟ᶜⁿ, Dd, model, p 
end
    
function invertSpectrum(𝐀, R, gf, model, CN, Dd, method)
	lb, ub = zeros(length(R)), ones(length(R))
	if hardbound == true
		ii = gf .< 1.0
		ub[ii] .= 0.0
	end	
 	x₀ = Normalize(R)
    
    Ax = Any[]	
	xλ = try 
		@match method begin
			"L₀DₓB" => @> invert(𝐀, R, LₖDₓB(0, 0.001, lb, ub)) 
            "LSQ₁"  => begin
                init = [1.0, 1.2]
				f(_, p) = (model(CN, p[1], Dd, p[2])).N
				fit = curve_fit(f, R, R, init, lower = zeros(2))
				Ax = fit.param
				kk = argmin(abs.(Ax[2] .- gf))
				@> zeros(length(gf)) setindex!(Ax[1], kk)
            end
            "LSQ₂"  => begin
                init = [0.5, 0.5, 1.1, 1.3]
                f(_, p) = (model(CN, p[1:2], Dd, p[3:4])).N
                fit = curve_fit(f, R, R, init, lower = zeros(4))
                Ax = fit.param
                kk = argmin(abs.(Ax[3] .- gf))
                ll = argmin(abs.(Ax[4] .- gf))
            @> zeros(length(gf)) setindex!(Ax[1], kk) setindex!(Ax[2], ll)
        end
			_ => throw("Not supported")
		end
	catch
		return lb
	end
end

function getInversion(𝐀, gf, ge, R, model, CN, Dd, i)
    xλ1 = invertSpectrum(𝐀, R, gf, model, CN, Dd, "L₀DₓB")
    xλ2 = invertSpectrum(𝐀, R, gf, model, CN, Dd, "LSQ₁")
    xλ3 = invertSpectrum(𝐀, R, gf, model, CN, Dd, "LSQ₂")

    r1 = @> RMSE(𝐀*xλ1, R) round(digits = 3)
    r2 = @> RMSE(𝐀*xλ2, R) round(digits = 3)
    r3 = @> RMSE(𝐀*xλ3, R) round(digits = 3)

    df1 = DataFrame(gf = htdma.gf[:,i], N = htdma.N[:,i], Color = "Measured")
    df3 = DataFrame(gf = gf, N = 𝐀*xλ1, Color = "𝐀<sub>2</sub>*L<sub>0</sub>D<sub>ε</sub>B ($(r1))")
    df4 = DataFrame(gf = gf, N = 𝐀*xλ2, Color = "𝐀<sub>2</sub>*LSQ<sub>1</sub> ($(r2))")
    df5 = DataFrame(gf = gf, N = 𝐀*xλ3, Color = "𝐀<sub>2</sub>*LSQ<sub>2</sub> ($(r3))")
    dfl = [df1;df3;df4;df5]

    dfr1 = DataFrame(gf = gf, Frequency = Normalize(xλ1), Color = "L<sub>0</sub>D<sub>ε</sub> B")
    dfr2 = DataFrame(gf = gf, Frequency = Normalize(xλ2), Color = "LSQ<sub>1</sub>")
    dfr3 = DataFrame(gf = gf, Frequency = Normalize(xλ3), Color = "LSQ<sub>2</sub>")
    dfr = [dfr1;dfr2;dfr3]

    p = plot_gf_dual(dfl, dfr)

    (method, xλ) = @match argmin([r1,r2,r3]) begin
        1 => ("L₀DₓB", xλ1)
        2 => ("LSQ₁", xλ2)
        3 => ("LSQ₂", xλ3)
        _ => "Error"
    end

    current = DataFrame(
        t   = htdma.timestamp[i], 
        Dd  = htdma.Dd[i],
        gl  = ge[2:end],
        gf  = gf, 
        gu  = ge[1:end-1],
        pdf = Normalize(xλ), 
        Method = method,
    )
    current, p
end

function runfile(htdmafile, smpsfile)
    global htdma = load_day_htdma(htdmafile)
    global smps = load_day_smps(smpsfile)

    n = length(htdma.timestamp)
    df = mapreduce(vcat, 1:n) do i
        𝐀, gf, ge, R, CN, Dd, model, p = prepInversion(i);
        current, p2 = getInversion(𝐀, gf, ge, R, model, CN, Dd, i);
        println(htdma.timestamp[i]," ", i,"/",n)
        current
    end    
end

const lpm = 1.6666666e-5

Λ,  δ  = getDMA(298, 950e2, 1lpm, 5lpm, 13.0, "TSI Long")
Λ₁, δ₁ = getDMA(298, 950e2, 0.63lpm, 5lpm, 0.0, "Brechtel")
Λ₂     = getDMA2(298, 950e2, 1lpm, 5lpm, 0.0, "Brechtel")

bins = 30
hardbound = true

# Run through several days. This is a bit of a hack and should be automated. 
# Since it's only a few days the effort wasn't worth it.
# Problem is there are sometimes multiple HTDMA files per day. 
hpath = "../../data/raw/sgpaoshtdma/"
spath = "../../data/raw/sgpaossmps/"
htdmafile = hpath*"sgpaoshtdmaE13.a1.20200206.000126.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200206.000000.nc"
dfa = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200207.000136.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200207.000000.nc"
dfb = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200208.000147.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200208.000000.nc"
dfc = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200209.000154.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200209.000000.nc"
dfd = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200210.000201.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200210.000000.nc"
dfe = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200211.000209.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200211.000000.nc"
dff = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200212.000216.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200212.000000.nc"
dfg = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200213.000223.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200213.000000.nc"
dfh = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200214.000231.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200214.000000.nc"
dfi = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200215.000240.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200215.000000.nc"
dfj = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200216.000249.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200216.000000.nc"
dfk = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200217.000258.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200217.000000.nc"
dfl = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200218.000307.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200218.000000.nc"
dfm = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200219.000316.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200219.000000.nc"
dfn = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200220.000325.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200220.000000.nc"
dfo = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200221.000334.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200221.000000.nc"
dfp = runfile(htdmafile, smpsfile)

htdmafile = hpath*"sgpaoshtdmaE13.a1.20200221.130649.custom.cdf"
smpsfile = spath*"sgpaossmpsE13.a1.20200221.000000.nc"
dfq = runfile(htdmafile, smpsfile)

[dfa; dfb; dfc; dfd; dfe; dff; dfg; dfh; dfi; dfj; dfk; dfl; dfm; dfn; dfo; dfp; dfq] |> 
    CSV.write("../../data/processed/sgpaoshtdmainverted.csv")
