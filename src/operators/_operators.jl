function kernel_from_integral(g, t)
    return g(t[2:end]) - g(t[1:end - 1])
end

function kernel_c2d(g, t)
    dt = similar(t)
    dt[1:end - 1] = diff(t)
    dt[end] = dt[end - 1]
    return g(t) .* dt
end # function

function response(g, u)
    y = similar(u)
    for t = 1:length(u)
        val = zero(y[t])
        for τ = 1:t - 1
            val += g[t - τ] * u[τ]
        end # for
        y[t] = val
    end # for
    return y
end # function

export response, kernel_c2d