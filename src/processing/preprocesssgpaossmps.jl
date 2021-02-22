# Preprocess SGP SMPS data. Reads data from all netcdf files in the directory
# and creates a common CSV format for storing size distribution data 

using CSV
using DataFrames

include("sgpLoaders.jl")

reshape_sgpsmps("../../data/raw/sgpaossmps/") |> CSV.write("../../data/processed/sgpaossmps.csv")
