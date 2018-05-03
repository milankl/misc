# LLVM example

f(k) = 5k - 1

function g(k)
    for i=1:10
        k = f(k)
    end
    return k
end

r = g(1)

code_native(g,(Int,))
