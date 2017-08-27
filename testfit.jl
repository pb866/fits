using LsqFit
using DataFrames

# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
model(x, p) = p[1]*exp.(-x.*p[2])

# some example data
# xdata: independent variables
# ydata: dependent variable
df = DataFrame()
df[:xdata] = linspace(0,10,20)
df[:ydata] = model(df[:xdata], [1.0 2.0]) + 0.01*randn(length(df[:xdata]))
p0 = [0.5, 0.5]
xdata = df[:xdata]
convert(Array{Number},xdata)
ydata = df[:ydata]
println(typeof(xdata))
fit = curve_fit(model, xdata, ydata, p0)
# fit is a composite type (LsqFitResult), with some interesting values:
#	fit.dof: degrees of freedom
#	fit.param: best fit parameters
#	fit.resid: residuals = vector of residuals
#	fit.jacobian: estimated Jacobian at solution

# We can estimate errors on the fit parameters,
# to get 95% confidence error bars:
errors = estimate_errors(fit, 0.95)

# The finite difference method is used above to approximate the Jacobian.
# Alternatively, a function which calculates it exactly can be supplied instead.
function jacobian_model(x,p)
    J = Array{Float64}(length(x),length(p))
    J[:,1] = exp.(-x.*p[2])    #dmodel/dp[1]
    J[:,2] = -x.*p[1].*J[:,1]  #dmodel/dp[2]
    J
end
fit = curve_fit(model, jacobian_model, xdata, ydata, p0)
println(fit)
