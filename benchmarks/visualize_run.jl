io = IOContext(stdout, :logbins => true)
for row in eachrow(results)
    show(io, MIME("text/plain"), row[2])
end