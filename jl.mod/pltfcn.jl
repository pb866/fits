module pltfcn

using DataFrames
using PyCall, PyPlot
import Winston
export plot_j

function plot_j(jvals,ifile,fit)

  # Initialise system time and output path/file name
  time = now()
  path = "./params_$(basename(splitext(ifile)[1]))"
  if isdir(path)
  else println("\033[93m'$path' created.\033[0m"); mkdir(path)
  end

  ptitle = beautify_chem(names(jvals))

  nj = length(fit)
  jplot = Array{PyCall.PyObject}(nj)
  xdata = jvals[1].*π./180.
  for i=1:nj
    p = floor(log10(fit[i].param[1]*exp(-fit[i].param[3])))
    l = fit[i].param[1]/10^p
    m = fit[i].param[2]
    n = fit[i].param[3]

    clf()
    jplot[i] = scatter(xdata,jvals[i+1]./(10^p),color=:black,marker = "s",label="TUV data")
    jplot[i] = plot(xdata,l.*(cos.(xdata)).^m.*exp.(-n./cos.(xdata)), color=:red, label="Parameterisation")
    size(2400, 1800)
    xticks([0.,0.1*π,0.2*π,0.3*π,0.4*π,0.5*π],
           [L"$0$",L"$\frac{1}{10}\pi$",L"$\frac{1}{5}\pi$",L"$\frac{3}{10}\pi$",L"$\frac{2}{5}\pi$",L"$\frac{1}{2}\pi$"])
    # xticks(0:π/8:π/2)
    # ymax(5.e-5)
    xlabel("χ")
    xlab = latexstring("j\\ /\\ 10^{$(convert(Int64,p))}\\,\\mathrm{s}^{-1}")
    ylabel(xlab)
    xlim(0.0,π/2)
    ylim(0.0,ceil(jvals[1,i+1]/10^p))
    # axis(ymin=0.0)
    grid(linestyle="dotted")
    # tight_layout()
    ltitle = latexstring(ptitle[i+1])
    title(ltitle)
    annotate("created $(Dates.format(time,"dd.mm.yyyy, HH:MM:SS"))",
    xy=[0.1;0.1],
    xycoords="axes fraction",
    ha="left",
    va="bottom")
    legend()
    savefig("$path/j$i.png", dpi=300)
  end

  # path = pwd()*"/py.fcn"
  # push!(pyimport("sys")["path"], path);
  # pyimport("testmod")[:multiplot](plot1,plot2)
  try
    run(`convert -quality 100 $path/'*'.png $path/plots.pdf`)
    cd(path)
    pfiles = filter(r"./*\.png",readdir())
    for i = 1:length(pfiles)  rm(pfiles[i],force=true)  end
    cd("..")
  catch
    println("\033[95mWARNING! UNIX tool 'convert' missing or unable to process plots.")
    println("No concatenation of plots in single pdf in results folder.")
    println("\033[96mMac users install imagemagick, e.g. using homebrew:")
    println("brew install imagemagick\033[0m")
  end

end #function plot_j

function beautify_chem(reactions)

  chem = String[]
  for rxn in reactions
    lstr = convert(String,rxn) # convert dataframe symbols to strings
    lstr = replace(lstr, "->", "\\stackrel{h\\nu}{\\longrightarrow}") # format arrows
    lstr = replace(lstr, " + hv", "") # remove + hv (now above formatted arrows)
    lstr = replace(lstr, r"([A-Za-z)])([0-9]+)", s"\1_{\2}") # make numbers in chemical formulas subscripts: {\1}_{\2}
    lstr = replace(lstr, r"\(([0-9]+)", s"(^{\1}") # except for atom's energetic state
    lstr = replace(lstr, ".", "^{.}") # raise radical dots
    lstr = replace(lstr, r"([0-9]+) ", s"\1\\,") # insert halfspace between stoichiometric indices and species
    lstr = replace(lstr, r"=", "\\!=\\!") # no spaces between double bonds
    lstr = replace(lstr, r"([A-Za-z0-9])-([A-Za-z])", s"\1\\!-\\!\2") # format dashes
    lstr = replace(lstr, r"C:", "{\\ddot C}") # format biradical functions
    lstr = replace(lstr, r"CH:", "{\\ddot C}H") # format biradical functions
    # define and correct exceptions:
    lstr = replace(lstr, "further products", "\\text{further products}") # transform to text
    # Ensure unslanted font
    lstr = "\\mathrm{"*lstr*"}"

    push!(chem, lstr)
  end

  return chem

end #function beautify_chem
end #module fitfcn
