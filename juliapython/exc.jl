# EXCEPTIONS

function squareroot(x)
    try
        return sqrt(x)
    catch
        return sqrt(x + 0im)
    end
end

squareroot(2)
squareroot(-2)
