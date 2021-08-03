# Script to evaluate the different inversion methods which are plotted in Figure S1
# This script is slow (1 hr runtime)

using DifferentialMobilityAnalyzers
using Distributions
using DataFrames
using MLStyle
using Lazy
using Random
using RegularizationTools
using ProgressMeter
using DataFrames
using CSV

include("../commonTDMAfunctions.jl")

""" function eval_error(x, Method, Phantom, Seed, Nt, Dd, xλ)
    Return the RMSE as a DataFrame
"""
function eval_error(x, Method, Phantom, Seed, Nt, Dd, xλ)
	RMSE = sqrt(sum((xλ .- x).^2.0)./length(x))
	return DataFrame(
		RMSE = RMSE, 
		Method = Method, 
		Phantom = Phantom, 
        Seed = Seed,
        Nt = Nt,
        Dd = Dd,
        bins = length(xλ)
    )
end

""" test_inv(seed, bins, Nt, Ddnm, TestCase)
    Simulates the the inversion for a given TestCase
"""
function test_inv(seed, bins, Nt, Ddnm, TestCase)
    Qcpc = 1.0 
    Dd = Ddnm*1e-9   

    Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, bins)
    Ax = [[0.66*Nt, 50.0, 1.4], [0.33*Nt, 130.0, 1.6]]
    𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
    gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, bins)
    model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 2.5, bins));

    edf = DataFrame[]
	dg = ge[1:end-1] .- ge[2:end]
    f = test_cases(TestCase, gf, ge, bins)
	println(sum(f))
	println(TestCase)
	println(Dd)
	println(bins)
    R = poisson_noise(1.0, 𝐀*f; seed = seed)
    x₀ = clean(R./sum(R))
    x₀[x₀ .>= 1] .= 0.999
    x₀[x₀ .<= 0] .= 0.001
    lb, ub = zeros(bins), ones(bins)

    map(0:2) do k
        @>> begin
  :q          invert(𝐀, R, LₖB(k, lb, ub))qq
            eval_error(f, "L$(k)B", TestCase, seed, Nt, Ddnm) 
            push!(edf)
        end
	end

	println("A")

    map(0:2) do k
        @>> begin
            invert(𝐀, R, LₖDₓB(k, 0.001, lb, ub)) 
            eval_error(f, "L$(k)DₓB", TestCase, seed, Nt, Ddnm) 
            push!(edf)
        end
	end
    
	println("B")

	map(0:2) do k
        @>> begin
            invert(𝐀, R, Lₖx₀B(k, x₀, lb, ub)) 
            eval_error(f, "L$(k)x₀B", TestCase, seed, Nt, Ddnm) 
            push!(edf)
        end
    end

	println("C")
    
	map(0:2) do k
        @>> begin
            invert(𝐀, R, Lₖx₀DₓB(k, x₀, 0.001, lb, ub)) 
            eval_error(f, "L$(k)x₀DₓB", TestCase, seed, Nt, Ddnm) 
            push!(edf)
        end
    end

	vcat(edf...)
end

# Conditions for the cases
seeds = collect(1000:1000:10000)
bins = collect(20:10:60)
Dd = [30, 50, 100, 200, 300]
Nt = [500, 1000, 5000, 50000]
phantom = ["Single Channel", "Two Channel", "Uniform", "Bimodal", "Truncated Normal"]

df = let 
    a = [DataFrame(s=i, bins=j, Dd=k, Nt=l, p=m) for i in seeds, j in bins, k in Dd, l in Nt, m in phantom] 
    vcat(a...)
end

mdf = @showprogress map(1:length(df[!,1])) do i
	a = try 
	    test_inv(df[i,:s], df[i,:bins], df[i,:Nt], df[i,:Dd], df[i,:p])
	catch
		nothing
	end
end

gdf = filter(x -> .~isnothing(x), mdf)

outdf = vcat(gdf...) |> CSV.write("../../data/processed/methodssummary.csv")
