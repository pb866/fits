photMCM
=======

Julia script to derive MCM parameterisations for photolysis processes from TUV output file
(in version 5.2 format).

    j / s-1 = l·(cos(χ))^m·exp(-n·sec(χ))


Dependencies
------------

The script is written in Julia (julialang). Compilers for all platforms and help
can be obtained at https://julialang.org/.

For python users, a python2.7 script is given as well. Python scripts are no longer
supported and run unreliably.


Running the julia script
------------------------

Place output file from TUV model run in the main folder and call julia script with

    julia photMCM.jl <scenario name>

The scenario name is the name defined by the parameter `outfil` in the TUV input file,
i.e. the file name without the `.txt` file ending.


Results
-------

The script creates a folder called `params_<scenario name>` with a text file
`parameters.dat` holding the parameters and 95% confidence intervals as well as
RMSE and the correlation factor R<sup>2</sup> for each reaction in the TUV output.

Each fit function is compared to the TUV calculated data in plots compiled in plots.pdf.


Differences in the python script
--------------------------------

In the python script specify the full file name with file ending of the TUV output file
as script argument.  
Output is given in the main folder as `<scenario name>.dat/pdf`.  
Another file `<scenario name>.par` is generated with the reaction numbers from the TUV
output file and the plain parameters for easier post-processing by further scripts.



Version history
===============

Version 2.0
-----------
- Script re-written in julialang
- Use file name without file ending as script input
- Output saved to separate folder `params_<scenario name>`
- Parameter file `\*.par` depricated


depricated – no longer supperted:
---------------------------------
Version 1.1
-----------
- Python2.7 script to derive parameters `l, m, n` for MCM photolysis processes
  from TUV (5.2) input file
- Text output with parameters and statistical information in <scenario name>.dat
  and plain parameters in <scenario name>.par
- Comparisons of parameterisations against TUV calculated data in plots in
  \<scenario name\>.pdf
