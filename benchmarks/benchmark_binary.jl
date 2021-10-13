# Julia script for benchmarking binary signal operators

using BenchmarkTools
using ProgressBars
using DataFrames
using JLD2

import Romeo

# Benchmark Convolution
# ---------------------
# func = Romeo.LTI.Operators.convolve
# filename = "benchmarks/runs/conv_N1000.jld2"

# Benchmark Convolution
# ---------------------
func = Romeo.LTI.Operators.deconvolve
filename = "benchmarks/runs/deconv_N1000.jld2"

types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64]
results = DataFrame(type=Type[], trial=BenchmarkTools.Trial[])
N = 1000
types_iter = ProgressBar(types)
for T in types_iter
    set_postfix(types_iter, type="$T")
    if T <: Integer
        g = rand(-10 * one(T):10 * one(T), N)
        u = rand(-10 * one(T):10 * one(T), N)
    else
        g = rand(T, N)
        u = rand(T, N)
    end
    r = @benchmark func($g, $u)
    push!(results, (T, copy(r)))
end

show(results)

jldsave(filename; results)