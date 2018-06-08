# pandigital number
p = [1,2,3,4,5,6,7,8,9]

function next_pandigital(p::Array{Int,1})
    N = 9
    i = N-1
    while i > 0
        # find trigger p[end-i]
        trigger = p[i]
        if trigger < p[i+1]
            s = p[i:end]    # subarray of numbers that have to be swapped
            a = minimum(s[s .> trigger])   # smallest number that is larger than p[end-i]
            sort!(filter!(x->x≠a,s))

            # do the swapping
            p[i] = a
            p[i+1:end] = s
            break
        else
            i -= 1
        end

        if i == 0   # pandigital number is 987654321
            p = p[end:-1:1]
        end
    end
    return p
end
