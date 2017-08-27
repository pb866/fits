#!/usr/local/bin/julia

push!(LOAD_PATH,"$(pwd())/jl.mod")

using rddata_df
using DataFrames

push!(ARGS,"")
jvals =  read_j(ARGS[1])
println(jvals)
