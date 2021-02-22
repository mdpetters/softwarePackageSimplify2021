using NetCDF
using DifferentialMobilityAnalyzers
using Dates
using Underscores
using MLStyle
using Lazy
using Interpolations
using DataFrames
using Gadfly
using Colors
using Printf

import Base.|>

|>(args...) = args[end](args[1:end-1]...)

struct SGP_HTDMA
    timestamp::Array{DateTime,1}
    dNdlog10D::Array{Float64,2}
    N::Array{Float64,2}
    Nt::Array{Float64,1}
    De1::Array{Float64,2}
    De2::Array{Float64,2}
    Dp::Array{Float64,2}
    ΔD::Array{Float64,2}
    Δlog10D::Array{Float64,2}
    Dd::Array{Float64,1}
    gf::Array{Float64,2}
    RH::Array{Float64,1}
    Q::Array{Float64,1}
end

struct SGP_SMPS
    timestamp::Array{DateTime,1}
    𝕟::Array{SizeDistribution,1}
    status_flag::Array{Int32,1}
    δ::DifferentialMobilityAnalyzer
    Λ::DMAconfig
end

function interpolateDataOntoδ(Dp, S, δ)
    itp = interpolate((Dp,), S, Gridded(Linear()))
    etp = extrapolate(itp, 0.0) 
    S = etp(δ.Dp)
    
    return SizeDistribution([],δ.De,δ.Dp,δ.ΔlnD,S./log(10),S.*δ.ΔlnD./log(10),:interpolated)
end

function load_day_smps(file)
    lpm = 1.6666666e-5
    Λ, δ = getDMA(298, 950e2, 1lpm, 5lpm, 13.0, "TSI Long")

    basetime = ncread(file, "base_time")
    bt = DateTime(1970, 1, 1, 0, 0, 0) + Dates.Second(convert(Int, basetime[1]))
    offset = ncread(file, "time_offset")
    timestamp = bt .+ Dates.Second.(offset)
    status_flag = ncread(file, "status_flag")
    psd = ncread(file, "number_size_distribution")
    Dp = ncread(file, "diameter_midpoint")
    Nt = ncread(file, "total_concentration")
    Qsh = ncread(file, "sheath_flow")
    psd[psd .== -9999] .= 0.0
    𝕟 = @_ map(interpolateDataOntoδ(Dp, psd[:,_], δ), 1:length(timestamp))

    if Qsh[1] ≠ 5.0
        throw("Error: wrong sheath flow")
    end

    SGP_SMPS(timestamp, 𝕟, status_flag, δ, Λ)
end

function load_day_htdma(file)
    basetime = ncread(file, "base_time")
    bt = DateTime(1970, 1, 1, 0, 0, 0) + Dates.Second(convert(Int, basetime[1]))
    offset = ncread(file, "time_offset")
    timestamp = bt .+ Dates.Second.(offset)
    dNdlog10D = ncread(file, "aerosol_concentration")
    Dp = ncread(file, "bin_center")
    ΔD = ncread(file, "bin_width")
    Dd = ncread(file, "dry_diameter_setting")
    RH = ncread(file, "humid_rh")
    Nt = ncread(file, "total_concentration")
    Q = ncread(file, "sample_flow")
    gf = hcat((@_ map(Dp[:, _] ./ Dd[_], 1:length(Dd)))...)

    De1 = Dp .- ΔD / 2
    De2 = Dp .+ ΔD / 2
    Δlog10D = log10.(De2 ./ De1)
    N = dNdlog10D .* Δlog10D

    SGP_HTDMA(timestamp, dNdlog10D, N, Nt, De1, De2, Dp, ΔD, Δlog10D, Dd, gf, RH, Q)
end

function sgp_to_csv(file, list; bins = false)
    ts = list.timestamp
    tint = Dates.value.(ts)
    𝕟 = list.𝕟[1]
    Dmin = minimum(𝕟.Dp)
    Dmax = maximum(𝕟.Dp)

    theDp = 𝕟.Dp
    upDp = 𝕟.De[1:end-1]
    lowDp = 𝕟.De[2:end]

    function single_line(j)
        df = DataFrame(File = file, timeISO8601 = ts[j], timeInt64 = tint[j])
        theN = @_ map(_ < 0 ? 0 : _, list.𝕟[j].N)
        map(i -> df[!, Symbol("DMA$i")] = [theN[i]], 1:length(𝕟.Dp))
        return df
    end
    df = map(single_line, 1:length(ts))... |> vcat

    dummy1 = DataFrame(
        File = "Lower bin bounds",
        timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
        timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
    )
    map(i -> dummy1[!, Symbol("DMA$i")] = [lowDp[i]], 1:length(𝕟.Dp))

    dummy2 = DataFrame(
        File = "Midpoints",
        timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
        timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
    )
    map(i -> dummy2[!, Symbol("DMA$i")] = [theDp[i]], 1:length(𝕟.Dp))

    dummy3 = DataFrame(
        File = "Upper bin bounds",
        timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
        timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
    )
    map(i -> dummy3[!, Symbol("DMA$i")] = [upDp[i]], 1:length(𝕟.Dp))


    if bins == true
        return [dummy1; dummy2; dummy3; df]
    else
        return df
    end
end

function getDMAdimensions(DMAtype)
    (r₁, r₂, l, geom) = @match DMAtype begin
        "TSI Long"      => (9.37e-3, 1.961e-2, 0.44369, :cylindrical)
        "High Flow DMA" => (0.05, 0.058, 0.6, :cylindrical)
        "RDMA"          => (2.4e-3, 50.4e-3, 10e-3, :radial)
        "Brechtel"      => (6.24e-2, 7.23e-2, 0.339, :cylindrical)
        _               => throw("Error")
    end
end

function getDMA(t, p, qsa, qsh, leff, DMAtype)
    (r₁, r₂, l, geom) = getDMAdimensions(DMAtype)
    Λ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, leff, :-, 6, geom)
    (Λ, setupDMA(Λ, dtoz(Λ, 600e-9), dtoz(Λ, 10e-9), 120))
end

function getDMA2(t, p, qsa, qsh, leff, DMAtype)
    (r₁, r₂, l, geom) = getDMAdimensions(DMAtype)
    Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, leff, :-, 6, geom)
end


function reshape_sgpsmps(path)
    files = readdir(path)
    list = load_day_smps(path * files[1])
    df1 = sgp_to_csv(files[1], list; bins = true)

	df2 = map(files[2:end]) do file
        @>> load_day_smps(path * file) sgp_to_csv(file)
    end

    return [df1; vcat(df2...)]
end

function df_to_stack(mdf, Dlow, Dup)
    df = mdf[4:end, :]
    Dp = convert(Vector, mdf[2, 4:end])
    ii = [BitArray([0, 0, 0]); (Dp .>= Dlow) .& (Dp .< Dup)]
    Dl = convert(Vector, mdf[1, ii])
    Dp = convert(Vector, mdf[2, ii])
    Du = convert(Vector, mdf[3, ii])

    n = length(Dp)
    mapreduce(vcat, 1:length(df[!, 1])-1) do i
        N = (convert(Vector, df[i, ii]))
        DataFrame(
            x1 = df[i, :timeInt64],
            x2 = df[i+1, :timeInt64],
            y1 = [Du[j] for j = 1:n],
            y2 = [Dl[j] for j = 1:n],
            N = N ./ maximum(N),
        )
    end
end
