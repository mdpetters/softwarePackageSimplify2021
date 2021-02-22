using Cairo
using Fontconfig
using Compose
using Gadfly
using NumericIO
using ColorTypes
using Printf
using CSV
using Dates
using Lazy
using MLStyle
using DataFrames

include("../commonTDMAfunctions.jl")

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

ylabel = log10.([10, 20, 50, 100, 200, 500])
lfuny = x -> ifelse(sum(x .== ylabel) == 1, @sprintf("%2i", exp10.(x)), "")
gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]

function plot_strip(ss, ee)
    label = ss:Dates.Day(1):ee
    xl0 = (Dates.value.(label), Dates.format.(label, "mm/dd"))
    xl1 = (Dates.value.(label), Dates.format.(label, ""))
    xl2 = map(i -> a = (i % 7 == true) ? (xl0[2])[i] : (xl1[2])[i], 1:length(label))
    xl = (Dates.value.(label), xl2)
    labeldict = Dict((xl[1])[i] => (xl[2])[i] for i = 1:length(label))

    Gadfly.plot(
        xmin = D[!, :x1],
        ymin = D[!, :y1],
        xmax = D[!, :x2],
        ymax = D[!, :y2],
        color = D[!, :N],
        Geom.rect,
        Guide.xlabel(""),
        Guide.ylabel("D (nm)", orientation = :vertical),
        style(
            plot_padding = [7.5mm, 5mm, 1mm, 0mm],
            key_position = :none,
            major_label_font = "PT Sans",
            minor_label_font = "PT Sans",
        ),
        Guide.xticks(ticks = Dates.value.(label)),
        Guide.yticks(ticks = log10.([10.0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500])),
        Scale.color_continuous(minvalue = 0, maxvalue = 0.99),
        Scale.y_log10(labels = lfuny),
        Scale.x_continuous(labels = i -> get(labeldict, i, "")),
        Coord.cartesian(xmin = Dates.value(ss), xmax = Dates.value(ee), ymin = log10(10.0)),
    )
end


df0 = CSV.read("../../data/processed/bbmlsmpsL0x0B.csv", DataFrame)
df2 = CSV.read("../../data/processed/bbmlsmpsL2B.csv", DataFrame)

D = df_to_stack(df0, 10, 600)
ss, ee = DateTime(2015, 1, 16, 0, 0, 0), DateTime(2015, 3, 6, 23, 59, 0)
p1 = plot_strip(ss, ee)

D = df_to_stack(df2, 10, 600)
p2 = plot_strip(ss, ee)

f1 = lch_diverge2()
xp = collect(0.5:0.002:0.7) .* w
n = length(xp)
yp = ones(n) .* 0.9h
colv = f1.(range(0, stop = 1.0, length = n))
dx = [xp[2] - xp[1] for i = 1:n]
dy = [0.08h for i = 1:n]

img1 = compose(
    context(),
    Compose.text(0.4w, 1h, "Normalized\ndN/dlnD (-)"),
    font("PT Sans"),
    Compose.text(0.50w, 0.85h, "0"),
    Compose.text(0.598w, 0.85h, "0.5"),
    Compose.text(0.698w, 0.85h, "1"),
    fill("grey30"),
    fontsize(8pt),
    (context(), Compose.rectangle(xp, yp, dx, dy), fill(colv)),
)

img2 = compose(
    context(),
    Compose.text(0.35w, 0.2h, "Date at midnight UTC in 2015 (month/day)"),
    font("PT Sans"),
    fill("grey30"),
    fontsize(10pt),
)

p = vstack(img1, p1, p2, img2)
set_default_plot_size(7inch, 12cm)
draw(PNG("figures/f06.png", dpi = 300), p)