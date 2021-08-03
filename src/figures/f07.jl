using TimeSeries
using Printf
using DataFrames
using Lazy
using CSV
using Colors
using Gadfly
using Compose
using MLStyle

import Cairo, Fontconfig

include("../commonTDMAfunctions.jl")

function get_plot_psd(ss, ee, D, Dl, Dh, yticks, ylabel, lab, pad)
    lfuny = x -> ifelse(sum(x .== ylabel) == 1, @sprintf("%2i", exp10.(x)), "")
    gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]

    label = ss:Dates.Hour(12):ee
    xl0 = (Dates.value.(label), Dates.format.(label, "m/d"))
    xl1 = (Dates.value.(label), Dates.format.(label, ""))
    xl2 = map(i -> a = (i % 2 == true) ? (xl0[2])[i] : (xl1[2])[i], 1:length(label))
    xl = (Dates.value.(label), xl2)
    labeldict = Dict((xl[1])[i] => (xl[2])[i] for i = 1:length(label))

    ii = .~isnan.(D[!, :N])
    p4 = plot(
        xmin = D[ii, :x1],
        ymin = D[ii, :y1],
        xmax = D[ii, :x2],
        ymax = D[ii, :y2],
        color = D[ii, :N],
        Geom.rect,
        Guide.xlabel(""),
        Guide.ylabel("D (nm)", orientation = :vertical),
        style(
            plot_padding = pad,
            major_label_font = "PT Sans",
            minor_label_font = "PT Sans",
            key_title_font_size = 8pt,
            key_position = :none
        ),
        Guide.xticks(ticks = Dates.value.(label)),
        Guide.yticks(ticks = yticks),
        Guide.colorkey(title = "$(lab)"),
        Scale.color_continuous(minvalue = 0, maxvalue = 1),
        Scale.y_log10(labels = lfuny),
        Scale.x_continuous(
            minvalue = Dates.value(ss),
            maxvalue = Dates.value(ee),
            labels = i -> get(labeldict, i, ""),
        ),
        Coord.cartesian(
            xmin = Dates.value(ss),
            xmax = Dates.value(ee),
            ymin = log10(Dl),
            ymax = log10(Dh),
        ),
    )
end

function get_plot_tdma(ss, ee, D, lab, pad)
    lfuny = x -> ifelse(sum(x .== ylabel) == 1, @sprintf("%2i", exp10.(x)), "")
    gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]

    label = ss:Dates.Hour(12):ee
    xl0 = (Dates.value.(label), Dates.format.(label, "m/d"))
    xl1 = (Dates.value.(label), Dates.format.(label, ""))
    xl2 = map(i -> a = (i % 2 == true) ? (xl0[2])[i] : (xl1[2])[i], 1:length(label))
    xl = (Dates.value.(label), xl2)
    labeldict = Dict((xl[1])[i] => (xl[2])[i] for i = 1:length(label))

    ii = .~isnan.(D[!, :pdf])
    p4 = plot(
        xmin = D[ii, :xmin],
        ymin = D[ii, :gl],
        xmax = D[ii, :xmax],
        ymax = D[ii, :gu],
        color = D[ii, :pdf],
        Geom.rect,
        Guide.xlabel(""),
        Guide.ylabel("gf <sub> $(lab) </sub> (-)", orientation = :vertical),
        style(
            plot_padding = pad,
            major_label_font = "PT Sans",
            minor_label_font = "PT Sans",
            key_title_font_size = 8pt,
            key_position = :none
        ),
        Guide.xticks(ticks = Dates.value.(label)),
        Guide.yticks(ticks = 1:0.2:1.6),
        Guide.colorkey(title = "$(lab)"),
        Scale.color_continuous(minvalue = 0, maxvalue = 1),
        Scale.x_continuous(
            minvalue = Dates.value(ss),
            maxvalue = Dates.value(ee),
            labels = i -> get(labeldict, i, ""),
        ),
        Coord.cartesian(
            xmin = Dates.value(ss),
            xmax = Dates.value(ee),
            ymin = 1,
            ymax = 1.6,
        ),
    )
end

function getDd(Dd)
    mdf = filter(:Dd => D -> D .== Dd, htdma_df)
    mdf = filter(:gf => g -> g .> 1, mdf)
    #mdf = filter(:Method => m -> m .== "LSQ₂", mdf)
    tuq = mdf[!,:t] |> unique |> sort
    mdf[!,:xmin] = Dates.value.(mdf[!,:t])
    mdf[!,:xmax] = Dates.value.(mdf[!,:t])

    map(1:length(tuq)-1) do i
        ii = mdf[!,:t] .== tuq[i]
        mdf[ii,:xmax] .= Dates.value(tuq[i+1])
        nothing
    end
    
    get_plot_tdma(ss, ee, mdf, "$(Dd)", [5mm,4.5mm,1mm,1mm])
end

sgp_df = CSV.read("../../data/processed/sgpaossmps.csv", DataFrame)
htdma_df = CSV.read("../../data/processed/sgpaoshtdmainverted.csv", DataFrame)
    
ss, ee = DateTime(2020, 2, 6, 0, 0, 0), DateTime(2020, 2, 22, 0, 0, 0)
subset(df) = [df[1:3, :]; filter(:timeISO8601 => t -> (t .> ss) .& (t .< ee), df)]

D1 = @> subset(sgp_df) df_to_stack(10, 550)

ylabel = log10.([10, 20, 50, 100, 200, 500])
yticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500])
p2 = get_plot_psd(ss, ee, D1, 10, 500, yticks, ylabel, "AOS SMPS", [5mm, 5mm, -10mm, 1mm])

p3 = getDd(50)
p4 = getDd(100)
p5 = getDd(150)
p6 = getDd(200)
p7 = getDd(250)

lch_diverge2 =
    function (l0 = 30, l1 = 100, c = 40, h0 = 260, h1 = 10, hmid = 20, power = 1.5)
        lspan = l1 - l0
        hspan1 = hmid - h0
        hspan0 = h1 - hmid
        function (r)
            r2 = 2r - 1
            return LCHab(
                min(80, l1 - lspan * abs(r2)^power),
                max(10, c * abs(r2)),
                (1 - r) * h0 + r * h1,
            )
        end
    end

f = lch_diverge2()
xp = collect(0.5:0.002:0.7) .* w
n = length(xp)
yp = ones(n) .* 0.4h
colv = f.(range(0, stop = 1.0, length = n))
dx = [xp[2] - xp[1] for i = 1:n]
dy = [0.08h for i = 1:n]

img1 = compose(
    context(),
    Compose.text(0.41w, 0.45h, "Colorscale"),
    font("PT Sans"),
    Compose.text(0.50w, 0.37h, "0"),
    Compose.text(0.588w, 0.37h, "0.5"),
    Compose.text(0.695w, 0.37h, "1"),
    fill("grey30"),
    fontsize(8pt),
    (context(), Compose.rectangle(xp, yp, dx, dy), fill(colv)),
)

img2 = compose(
    context(),
    Compose.text(0.38w, 0.2h, "Date at midnighy UTC (month/day)"),
    font("PT Sans"),
    fill("grey30"),
    fontsize(10pt),
)

set_default_plot_size(6.5inch, 6.5inch)
p = vstack(img1, p2, p7, p6, p5, p4, p3, img2)
Gadfly.draw(PNG("figures/f07.png", dpi = 300), p)
