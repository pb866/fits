"""
Module pltfcn
=============
version 1.1
-------------

List of function focused on plotting and writing file output.
contains:
    - scatdat
    - pfit
    - texstr
"""

import sys
#______________________________________________________________________________________________

def scatdat(rxn,data):
    """
        Function scatdat
        ================

        Purpose:
            Derive matrix with plot data and vector with associated photoreactions.

        Variables:
        I/O
            rxn:        matrix with indices and labels of available photoreactions
            data:       matrix with sza-dependent j values from TUV model
            om:         matrix with all orders of magnitudes (for sza column 0 dummy inserted)

        internal:
            h,x,y:      counters/indices
            o1,o2:      order of magnitude of maximum value in each j data column

        Dependencies:
        uses:           numpy, datfcn.order
        called from:    photMCM (main)
        """

    # initilise data
    import numpy as np
    from datfcn import order
    om    = []

    # determine order of magnitude from maximum (first) value for nicer output
    for y in range(len(rxn)):
        o1, o2 = order(data[0][y+1])
        om = np.append(om,o2)
    om = np.insert(om,0,1.)
        # write data matrix
    for x in range(len(data)):
        for y in range(len(rxn)+1):
            if (y == 0):
                data[x][0] = np.deg2rad(data[x][0])
    return rxn, data, om
#______________________________________________________________________________________________

def pfit(fout,fpar,pp,rxn,xdata,ydata,y,p,om,l,m,n,ol,el,em,en,rsquared,ss_tot,rmse,scen):
    """
    Function pfit
    =============

    Purpose:
        Produce output:
            - Data file with optimised parameters and statistical data
            - Matrix with parameters for further processing
            - Plots with TUV data and fitted functions

    Variables:
    I/O:
        fout,fpar:      identifiers for output files
        rxn:            matrix with indices and labels of available photoreactions
        p:              matrix with all optimised paramters from least square fit
        l,m,n:          optimised paramters from least square fit
        el,em,en:       confidence of parameters
        ol:             order of magnitude of parameter l
        om:             matrix with all orders of magnitudes (for sza column 0 dummy inserted)
        rsquared:       correlation coefficient
        ss_tot:         standard deviation
        rmse:           root mean square error
        y:              index (for addressing positions in rxn-matrix)
        scen:           scenario name

    internal:
        now:            variable for present time
        x:              x data for function plotting
        fig,axes,
        ymax,ymin:      variables used for plotting/formatting plots

    Dependencies:
        uses:           datetime, numpy, matplotlib.pyplot, texstr
        called from:    datfcn.xydat
    """

    import numpy as np
    from datetime import datetime as dt
    import matplotlib.pyplot as plt
    x = np.linspace(0.,np.pi/2,100) # define x-data for function plotting
    now = dt.now() # set timestamp

    # write data in table in file <scen>.dat
    fout.write("(%.3f%s%.3f)E%i\t %.3f%s%.3f \t %.3f%s%.3f \t  %.3e \t %.4f \t%s\n"
               %(l,u"\u00B1",el,ol,m,u"\u00B1",em,n,u"\u00B1",en,rmse,rsquared,rxn[y][1]))

    # write file with parameter matrix for further processing
    fpar.write("j%i\t%.3e\t%.3f\t%.3f\n" % (rxn[y][0],p[0],m,n,))

    # create plots with TUV data and fitted functions
    fig = plt.figure(y)
    fig.set_size_inches(6.,4.)
    axes = plt.gca()
    plt.xlabel('sza')
    tex = texstr(rxn[y][1])

    plt.ylabel('$j\mathrm{[%s]\ /\ 10^{%i}\,s^{-1}}$' % (tex, int(np.log10(om[y+1]))), fontsize=10)
    plt.xticks([0.,0.1*np.pi,0.2*np.pi,0.3*np.pi,0.4*np.pi,0.5*np.pi],
               [r'$0$',r'$\frac{1}{10}\pi$',r'$\frac{1}{5}\pi$',r'$\frac{3}{10}\pi$',r'$\frac{2}{5}\pi$',r'$\frac{1}{2}\pi$'])
    plt.plot(xdata,ydata/om[y+1],'ks')
    plt.plot(x,p[0]/om[y+1]*np.cos(x)**m*np.exp(-n/np.cos(x)),'r-', linewidth=2.0)
    ymin,ymax = axes.get_ylim()
    plt.text(0.03*np.pi,0.1*ymax,'created: %s.%s.%s, %s:%02d' % (now.day,now.month,now.year,now.hour,now.minute), fontsize=9)
    pp.savefig(fig)
    plt.close(fig)

    return None

#______________________________________________________________________________________________

def texstr(orig):
    """
        Function texstr
        ===============

        Purpose:
            Format string for nicer axis labels using LaTeX nomenclature.
                - Arrows used for -> (with hv above)
                - numbers as subscripts with following exceptions:
                    - O(1D)/O(3P)
                - biradical dots above C
                - no spacing between = of double bonds
                - further corrections of spacing

        Variables:
        I/O:
        orig:       original string
        repl:       formatted string

        Dependencies:
        uses:           re
        called from:    pltfcn.pfit
        """

    import re

    repl = orig.replace('->',r'\stackrel{h\nu}{\longrightarrow}') # format arrows
    repl = repl.replace(' + hv',r'') # remove + hv (now above formatted arrows)
    repl = re.sub('([A-Za-z)])([0-9]+)',r'\1_{\2}',repl) # make numbers in chemical formulas subscripts: {\1}_{\2}
    repl = re.sub('\(([0-9]+)',r'(^{\1}',repl) # except for atom's energetic state
    repl = re.sub('\.','^{.}',repl) # raise radical dots
    repl = re.sub('([0-9]+) ',r'\1\,',repl) # insert halfspace between stoichiometric indices and species
    repl = re.sub(r'=',r'\!=\!',repl) # no spaces between double bonds
    repl = re.sub(r'([A-Za-z0-9])-([A-Za-z])',r'\1\!-\!\2',repl) # format dashes
    repl = re.sub(r'C:',r'{\ddot C}',repl) # format biradical functions
    repl = re.sub(r'CH:',r'{\ddot C}H',repl) # format biradical functions
    # define and correct exceptions:
    repl = re.sub('further products','further\ products',repl) # include space

    return repl

#______________________________________________________________________________________________
