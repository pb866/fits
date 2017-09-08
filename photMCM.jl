#!/usr/local/bin/julia


"""
# Module photMCM

Derive fit parameters from TUV 5.2.1 output and list parameters and statistical data
in text file. Plot TUV and fit data.
"""
module photMCM

# Add path of self-made modules
push!(LOAD_PATH,"$(pwd())/jl.mod")

println("initialise...")
# Import Libraries/Modules
using DataFrames
using rddata, fhandle, fitfcn, pltfcn

# Initialise system time and output path/file name
time = now()

# Append ARGS by empty strings
# to avoid error messages in case of missing input
push!(ARGS,"")

# Read dataframe with j values from TUV output file
println("load data...")
ifile, iofolder = get_fnames(ARGS)
jvals =  read_j(ifile)
# Derive parameterisations for j values
println("fit data...")
fit, sigma, rmse, R2 = fit_j(jvals)
# Generate output
println("plot data...")
plot_j(jvals,iofolder,fit,time)
wrt_params(names(jvals),fit,sigma,rmse,R2,iofolder,time)

end #module photMCM
