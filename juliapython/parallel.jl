# PARALLEL example

a = SharedArray{Float64}(10)
@parallel for i = 1:10
    a[i] = i
end

println(a)
