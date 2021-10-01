using JLD2
using DataFrames
using BenchmarkTools

file_names = ["benchmarks/runs/conv_N1000.jld2", "benchmarks/runs/deconv_N1000.jld2"]

run = jldopen(file_names[1])
results_conv = run["results"]

run = jldopen(file_names[2])
results_deconv = run["results"]

results = DataFrame()
insertcols!(results, 1, :type => results_conv[!, :type], :conv => results_conv[!, :trial], :deconv => results_deconv[!, :trial])

# io = IOContext(stdout, :logbins => true)
# for row in eachrow(results)
#     show(io, MIME("text/plain"), row[2])
# end