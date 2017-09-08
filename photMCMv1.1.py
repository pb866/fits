#!/Users/peter/anaconda/bin/python

"""
####################
#                  #
#  photMCM (main)  #
#  version 1.1     #
#                  #
####################

Purpose:
    Script to read in TUV output and calculate parameterisations for photolysis reactions in MCMv4.0

Instructions:
    Run script with command:

        python photMCM <input file> <output>

    where input file is the name of the file from the TUV5.2.1 model runs and
    output is the name for the files with the optimized parameters (<output>.dat)
    and the plots of the TUV data and the fitted functions (<output>.ps).

    If left blank, output will have the same names as the input assuming the file
    ends with '.txt. Otherwise, output will be set to 'photMCM'. User input is required
    for blank input file during the application of the script.

    To make script run on a different maschine, copy script, the folder 'py.fcn/v1.1' with the modules
    datfcn, pltfcn, and fitfcn and change folder path of sys.path on l. 78 in the main script.

Input data:
    Use the output of TUV5.x (or model output with equivalent output) as provided by the TUV outfil.
    Calculations have to be performed sza rather than time dependent
    (Loop over solar zenith angles, set 'lzenit' true in TUV).


Variables:
    ifile:          identifier for input file with TUV photlysis data
                    (read from command line or from user input)
    spath           path of the modules with system user added
    fout,fpar:      file identifiers for output temporary files
    rxn:            matrix with available photoreactions and indices
    data:           matrix with sza-dependent j values from TUV
    om:             list of maximum order of magnitudes for l-parameters in MCM parameterisation
    scen:           scenario name for output files (from command line or photMCM)
    str:            string for rename ps file with plots to scenario name
    now:            variable for present time

Dependencies:
    uses:           sys,os,datetime,datfcn,pltfcn.scatdat

for help, see also:
    http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.curve_fit.html
    http://www.walkingrandomly.com/?p=5215

Author:

    Please report any bugs/issues/comments or direct your questions to:

    Peter Brauer
    Wolfson Atmospheric Chemistry Laboratories
    Department of Chemistry
    University of York
    York, YO10 5DD

    email: peter.brauer@york.ac.uk
    phone: +44 (0)1904 32 4758
    skype: p.braeuer


This script may be used, redistributed and/or altered for non-commercial purposes
under the GNU COMMON USER licence.
"""


# Load system functions
import sys
reload(sys)  # Reload does the trick!
sys.setdefaultencoding('UTF8') # output of UNICODE characters
import os
sys.path.insert(0,os.path.expanduser("~/Google Drive/Applications/data/photMCM/py.fcn/v1.1"))
from datetime import datetime as dt

#load own functions
from datfcn import *
from pltfcn import scatdat


# Load TUV data
try:
    ifile = sys.argv[1]
except:
    ifile = raw_input("Enter file with TUV data: ")
# retrival of fitting data
rxn, data = load(ifile)
data = trans(data)

# Initialise output
# text file with optimised parameters and statistical data:
try:
    scen = sys.argv[2]
except:
    if (ifile[-4:]=='.txt'):
        scen = ifile[:-4]
    else:
        scen = 'photMCM'

now = dt.now()
fout = open("%s.dat" % scen,'w+')
fout.write("Parameters and statistical data for parameterisation of photolysis processes in MCMv4.0.\n")
fout.write("j / s-1 = l%s(cos(x))^m%sexp(-n%ssec(x))\n\n" % (u"\u00B7",u"\u00B7",u"\u00B7"))
fout.write("created %s.%s.%s, %s:%s\n\n\n" % (now.day,now.month,now.year,now.hour,now.minute))
fout.write("               P a r a m e t e r s        \t\t\t  S t a t i s t i c s\n")
fout.write("    l / s-1     \t      m      \t      n      \t RMSE / s-1 \t   R^2  \tReaction\n")


# plot file for further gnuplot processing:
fpar = open("%s.par" % scen,'w+')
# write data for scatter plots in gnuplot input file
rxn, data, om = scatdat(rxn,data)

# info on screen
print "Done loading data.\nStart calculating parameterisations.\n\n"
print "Working on:\n"

# loop over photoreactions and curve fitting
xydat(data,rxn,90.,om,fout,fpar,scen)

# close all open files
fout.close()
fpar.close()


# info screen
print "\n\nDone.\n"
print "See file \'%s.dat\' for optimised parameters" % scen
print "of parmaterisation for photolysis in MCMv4.0 and statistical data."
print "Plots of TUV calculated data and least square fits are provided in \'%s.ps\'" % scen
print "Matrix with parameters is provided in \'%s.par\' for further processing." % scen



"""
plt.plot(xdata,ydata,xdata,phot(p,xdata))
plt.title(rxn[y][1])
plt.legend("TUV output","Fit")
os.system('mkdir plots')
fig = "plots/j"+str(rxn[y][0])+".eps"
plt.savefig(fig, format='eps', dpi=300)
"""
