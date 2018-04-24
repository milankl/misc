function box_count(x,y,z)
    #= BOX COUNT algorithm, assumes a long time series of x,y,z that are normalized in the (0,1) cube. =#

    # for n in nvec: n = 1/ϵ where n is the number of boxes in one direction, ϵ is the scale
    nvec = 2.^(3:8)
    Nmax = zeros(nvec)

    # largest distance between consequtive points - to compare with box size ϵ
    #Δx = maximum(sqrt.((x[1:end-1]-x[2:end]).^2 + (y[1:end-1]-y[2:end]).^2 + (z[1:end-1]-z[2:end]).^2))

    Nt = length(x)

    for (ni,n) in enumerate(nvec)
        # preallocate
        B = zeros(Int,n,n,n)    # could also be set up as boolean, indicates whether a box has been visited or not
        N = zeros(Int,Nt+1)     # counter array, set up as array to also check for convergence behaviour if necessary

        for i = 1:Nt
            xi = Int(ceil(x[i]*n))  # allocate a point into its box via rounding
            yi = Int(ceil(y[i]*n))
            zi = Int(ceil(z[i]*n))
            N[i+1] = N[i]+(1-B[xi,yi,zi])   # only add to counter if box was previously not visited
            B[xi,yi,zi] = 1     # mark box as visited
        end

        Nmax[ni] = N[end]
    end

    # estimate slope via linear regression
    interscept,slope = linreg(-log10.(1./nvec),log10.(Nmax))

    return slope
end
