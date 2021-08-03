using Gadfly
using CSV
using Statistics
using StatsBase
using DataFrames

df = CSV.read("../../data/processed/methodssummary.csv", DataFrame)

p = plot(
    df,
    xgroup = :Phantom,
    ygroup = :bins,
    Geom.subplot_grid(
        layer(x = :Method, y = log10.(df[!,:RMSE]), Geom.boxplot(;suppress_outliers=true), Theme(default_color = "steelblue3")),
    ),
    Guide.xlabel("Method by test case"),
    Guide.ylabel("Log10 root mean square error by number of bins"),
    Scale.y_continuous(minvalue=-3, maxvalue=0)
    )
img = SVG("figures/fs01.svg", 9inch, 6.5inch)

draw(img, p)
