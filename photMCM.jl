#!/usr/local/bin/julia

push!(LOAD_PATH,"$(pwd())/jl.mod")

println("initialise...")
# Import Libraries/Modules
using DataFrames
using rddata, fitfcn, pltfcn

# Append ARGS by empty strings
# to avoid error messages in case of missing input
push!(ARGS,"")

# Read dataframe with j values from TUV output file
println("load data...")
jvals =  read_j(ARGS[1])
# Derive parameterisations for j values
println("fit data...")
fit, sigma, rmse, R2 = fit_j(jvals)
# Generate output
println("plot data...")
plot_j(jvals,fit)
