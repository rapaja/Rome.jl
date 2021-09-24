using Romeo
using Plots

α = 0.5
β = 0.5
Δt = 0.01
N = 100

t = Δt:Δt:((N - 1) * Δt)
g_c = t -> Romeo.Fractional.ComplexOrder.kernel_fcn_co_integral(t, α, β)
g = Romeo.Operators.kernel_c2d(g_c, t)
g_i = t -> Romeo.Fractional.ComplexOrder.kernel_fcn_co_integral(t, α + 1, β)
gi = Romeo.Operators.kernel_from_integral(g_i, t)
u = ones(size(t))
y = Romeo.Operators.response(g, u)
yi = Romeo.Operators.response(gi, u)
ye = Romeo.Fractional.ComplexOrder.kernel_fcn_co_integral(t, α + 1, β)

resp = [y ye yi]

# Plots.plot(t, g)
Plots.plot(t, resp)