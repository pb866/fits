__precompile__()


"""
# Module _fhandle_

Routines for file handling, reading/storing input data and writing output data
to text files.

# Functions

- get_fnames (public)
- wrt_params (public)
"""
module fhandle
using DataFrames
export get_fnames, wrt_params


"""
    get_fnames(scen)

From the scenario name (scen) used in TUV derive names for the TUV output file
and the output folder.
"""
function get_fnames(scen)

  # Test script argument and ask for input file, if missing
  if scen == ""
    scen = print("Enter file with TUV scenario name")
    println("(Name of output file without '.txt'): ")
    scen = readline(STDIN)
  else
    scen = ARGS[1]
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  try mkdir(iofolder)
  catch
    print("\033[95mFolder '$iofolder' already exists! Overwrite ")
    print("(\033[4mY\033[0m\033[95mes/\033[4mN\033[0m\033[95mo)?\033[0m ")
    confirm = readline(STDIN)
    if lowercase(confirm) == "yes" || lowercase(confirm) == "y"
      cd(iofolder); files = readdir(); for file in files  rm(file)  end; cd("..")
    else println("Programme aborted. Exit code '98'."); exit(98)
    end
  end

  # Define TUV file
  ifile = scen*".txt"
  if isfile(ifile) == false
    println("$ifile does not exist. Script terminated. Exit code: 99."); exit(99)
  end

  # return file and folder names
  return ifile, iofolder
end #function get_fnames


"""
    wrt_params(rxn, fit, sigma, rmse, R2, iofolder, time)

For each reaction (rxn), from param in fit and sigma, print parameters
and 95% confidence together with RMSE and R^2 to file 'parameters.dat'
in the designated output folder (iofolder).

Print the time of creation (time) to the output file.
"""
function wrt_params(rxn, fit, sigma, rmse, R2, iofolder, time)

  #transform dataframe column symbols to strings
  rxn = convert.(String,rxn)

  # Open output file
  open("$iofolder/parameters.dat","w") do f
    # Print header
    println(f,
    "Parameters and statistical data for parameterisation of photolysis processes")
    println(f,
    "in the Master Chemical Mechanism (MCM; http://mcm.leeds.ac.uk/MCMv3.3.1/) using:")
    println(f, "\nj / s-1 = l·(cos(x))^m·exp(-n·sec(x))")
    println(f, "\n\ncreated $(Dates.format(time,"dd.mm.yyyy, HH:MM:SS"))")
    println(f,"\n                P a r a m e t e r s               S t a t i s t i c s")
    println(f,"    l / s-1              m              n         RMSE / s-1    R^2      Reaction")

    # Loop over reactions
    for i = 1:length(fit)
      # Determine order of magnitude of parameter l
      p = floor(log10(fit[i].param[1]))
      # Print parameters, statistical data, and reaction label to output file
      @printf(f,"(%.3f±%.3f)E%d    %.3f±%.3f    %.3f±%.3f    %.3e    %.4f    %s\n",
      fit[i].param[1]/10^p,sigma[i][1]/10^p,p,
      fit[i].param[2],sigma[i][2],fit[i].param[3],sigma[i][3],rmse[i],R2[i],rxn[i+1])
    end
  end
end # function wrt_params
end # module fhandle
