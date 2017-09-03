module pltfcn

using DataFrames
using PyPlot
import Winston
export plot_j

function plot_j(jvals,fit)

  PyPlot.svg(true)
  time = now()
  path = "./results"
  if isdir(path)
  else
    println("\033[93m'$path' created.\033[0m"); mkdir(path)
  end

  xdata = jvals[1].*π./180.
  l = fit[1].param[1]*1.e5
  m = fit[1].param[2]
  n = fit[1].param[3]

  plot1 = scatter(xdata,jvals[2]./1.e-5,color=:black,marker = "s",label="TUV data")
  plot1 = plot(xdata,l.*(cos.(xdata)).^m.*exp.(-n./cos.(xdata)), color=:red, label="Parameterisation")
  size(2400, 1800)
  xticks([0.,0.1*π,0.2*π,0.3*π,0.4*π,0.5*π],
         [L"$0$",L"$\frac{1}{10}\pi$",L"$\frac{1}{5}\pi$",L"$\frac{3}{10}\pi$",L"$\frac{2}{5}\pi$",L"$\frac{1}{2}\pi$"])
  # xticks(0:π/8:π/2)
  # ymax(5.e-5)
  xlabel("χ")
  ylabel(L"j / \mathrm{s}^{-1}")
  xlim(0.0,π/2)
  ylim(0.0,5.0)
  # axis(ymin=0.0)
  grid()
  tight_layout()
  annotate("created $(Dates.format(time,"dd.mm.yyyy, HH:MM:SS"))",
  xy=[0.1;0.1],
  xycoords="axes fraction",
  ha="left",
  va="bottom")
  legend()
  savefig("$path/test.png", dpi=300)

  clf()
  xdata = jvals[1].*π./180.
  l = fit[2].param[1]*1.e4
  m = fit[2].param[2]
  n = fit[2].param[3]

  plot2 = scatter(xdata,jvals[3]./1.e-4,color=:black,marker = "s",label="TUV data")
  plot2 = plot(xdata,l.*(cos.(xdata)).^m.*exp.(-n./cos.(xdata)), color=:red, label="Parameterisation")
  size(2400, 1800)
  xticks([0.,0.1*π,0.2*π,0.3*π,0.4*π,0.5*π],
         [L"$0$",L"$\frac{1}{10}\pi$",L"$\frac{1}{5}\pi$",L"$\frac{3}{10}\pi$",L"$\frac{2}{5}\pi$",L"$\frac{1}{2}\pi$"])
  # xticks(0:π/8:π/2)
  # ymax(5.e-5)
  xlabel("χ")
  ylabel(L"j / \mathrm{s}^{-1}")
  xlim(0.0,π/2)
  ylim(0.0,5.0)
  # axis(ymin=0.0)
  grid()
  tight_layout()
  annotate("created $(Dates.format(time,"dd.mm.yyyy, HH:MM:SS"))",
  xy=[0.1;0.1],
  xycoords="axes fraction",
  ha="left",
  va="bottom")
  legend()
  savefig("$path/test2.png", dpi=300)
  # path = pwd()*"/py.fcn"
  # push!(pyimport("sys")["path"], path);
  # pyimport("testmod")[:multiplot](plot1,plot2)
  try
    run(`convert -quality 100 $path/'*'.png $path/results.pdf`)
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
end #module fitfcn
