# PARALLEL example 2
N = 200000000
nheads = @parallel (+) for i = 1:N
    Int(rand(Bool))
end

println(nheads/N)
