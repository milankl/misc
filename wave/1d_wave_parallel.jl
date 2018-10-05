using PyPlot
using MPI

# OPERATORS - PERIODIC
function ∂x!(du::AbstractVector,u::AbstractVector)
    m = length(du)
    @boundscheck m+1 == length(u) || throw(BoundsError())

    @inbounds for i ∈ 1:m
        du[i] = one_over_dx*(u[i+1] - u[i])
    end
end

function Ix!(du::AbstractVector,u::AbstractVector)
    m = length(du)
    @boundscheck m+1 == length(u) || throw(BoundsError())

    @inbounds for i ∈ 1:m
        du[i] = 0.5*(u[i] + u[i+1])
    end
end


# Forcing
F(x,t) = -F0*sin.(8π*x/Lx .+ 200*t).*sin.(π*x/Lx).^2/rho

# RIGHT HAND SIDE
function rhs!(du,dη,u,η,h,u²,h_u,u_h,dudx,p,dpdx,U,dUdx,t)

    @views h .= η .+ H
    @views u² .= u.^2

    Ix!(h_u,h)
    Ix!(u_h,u²)
    ∂x!(dudx,u)

    # Bernoulli potential + stress "tensor"
    @views p .= -.5*u_h .- g*η[2:end] .+ ν*dudx
    #@views p .= -g*η[2:end] .+ ν*dudx

    # momentum
    ∂x!(dpdx,p)
    @views du[2:end-1] .= dpdx .+ F(x_up,t)./h_u[2:end]

    # continuity
    @views U .= u[1:end-1].*h_u
    ∂x!(dUdx,U)
    @views dη[2:end-1] .= -dUdx
end

# grid points per process
function subdomain_size(N::Int,size::Int,rank::Int)
    # rank is assumed to count up from 0
    if mod(N,size) == 0
        # domain size is divisible by the number of processors
        Np = N ÷ size
    else
        throw(error("Domain size must be divisible by number of processes."))
        # if rank+1 < size
        #     # the first n-1 processors get a slightly larger chunk of the domain
        #     Np = N ÷ size + 1
        # elseif rank + 1 == size
        #     # the last processor gets what's left
        #     Np = N - (size-1)*(N÷size + 1)
        # else
        #     throw(error("Rank higher than number of processors."))
        # end
    end
    return Np
end

function neighbours(size::Int,rank::Int)
    # left, right neighbour
    lnb,rnb = mod(rank-1,size),mod(rank+1,size)
end

function add_halo(u::AbstractVector,η::AbstractVector,rmsg::AbstractVector,lmsg::AbstractVector)
    u = cat(0.,u,0.,dims=1)
    η = cat(0.,η,0.,dims=1)

    ghost_points!(u,η,rmsg,lmsg)
    return u,η
end

function ghost_points!(u::AbstractVector,η::AbstractVector,rmsg::AbstractVector,lmsg::AbstractVector)
    if size == 1
        # sequential
        u[1] = u[end-1]
        u[end] = u[2]

        η[1] = η[end-1]
        η[end] = η[2]

    else
        # parallel with communication
        MPI.Isend([u[2],η[2]],lnb,rank+10,comm)
        MPI.Isend([u[end-1],η[end-1]],rnb,rank+20,comm)

        rreq = MPI.Irecv!(rmsg,rnb,rnb+10,comm)
        lreq = MPI.Irecv!(lmsg,lnb,lnb+20,comm)

        MPI.Waitall!([rreq,lreq])

        u[end],η[end] = rmsg
        u[1],η[1] = lmsg
    end
end

function preallocate(u::AbstractVector,η::AbstractVector)
    # u, η already have the halo

    u0,η0 = zero(u),zero(η)
    u1,η1 = zero(u),zero(η)
    du,dη = zero(u),zero(η)

    u²,h = zero(u),zero(η)
    p,dpdx = zeros(length(η)-1),zeros(length(η)-2)

    h_u = zeros(length(η)-1)
    u_h = zeros(length(u)-1)
    dudx = zeros(length(u)-1)

    U = zeros(length(u)-1)
    dUdx = zeros(length(u)-2)

    return u0,η0,u1,η1,du,dη,h,u²,h_u,u_h,dudx,p,dpdx,U,dUdx
end

function time_integration(Nt,u,η)

    # parallel: message from right, left
    rmsg = Array{Float64}(undef,2)
    lmsg = Array{Float64}(undef,2)

    # pre-allocate memory
    u,η = add_halo(u,η,rmsg,lmsg)
    u0,η0,u1,η1,du,dη,h,u²,h_u,u_h,dudx,p,dpdx,U,dUdx = preallocate(u,η)

    t = 0.

    # for output
    u_out = zeros(Nt+1,length(u)-2)     # -2 to remove ghost points
    η_out = zeros(Nt+1,length(η)-2)

    # store the initial conditions
    u_out[1,:] = u[2:end-1]
    η_out[1,:] = η[2:end-1]

    for i = 1:Nt

        ghost_points!(u,η,rmsg,lmsg)
        u1 .= u
        η1 .= η

        for RKi = 1:4
            if RKi > 1
                ghost_points!(u1,η1,rmsg,lmsg)
            end

            rhs!(du,dη,u1,η1,h,u²,h_u,u_h,dudx,p,dpdx,U,dUdx,t)

            if RKi < 4 # RHS update for the next RK-step
                u1 .= u .+ RKβ[RKi]*dt*du
                η1 .= η .+ RKβ[RKi]*dt*dη
            end

            # Summing all the RHS on the go
            u0 .+= RKα[RKi]*dt*du
            η0 .+= RKα[RKi]*dt*dη

        end

        u .= u0
        η .= η0
        t += dt

        # # output
        u_out[i+1,:] .= u[2:end-1]
        η_out[i+1,:] .= η[2:end-1]

    end
    u_out,η_out
end

MPI.Init()
comm = MPI.COMM_WORLD

const rank = MPI.Comm_rank(comm)
const size = MPI.Comm_size(comm)

const g = 10.
const H = 10.
const N = 500
const Nt = 3000
const cfl = 0.9
const rho = 1.
const ν = 0.006
const F0 = 500.
const Lx = 1.
const dx = Lx/N
const one_over_dx = 1/dx
const cph = sqrt(g*H)
const dt = cfl * dx / cph
const RKα = [1/6.,1/3.,1/3.,1/6.]
const RKβ = [0.5,0.5,1.]

# grid
const x_h = dx/2:dx:Lx
const x_u = dx:dx:Lx        # no u-point at x=0 but at x=L (periodicity)

# subdomains for parallel
const Np = subdomain_size(N,size,rank)
const x_up = x_u[(rank*Np+1):(rank+1)*Np]
const lnb,rnb = neighbours(size,rank)

# initial conditions
u_ini = fill(0.,Np)
η_ini = fill(0.,Np)

@time u,η = time_integration(Nt,u_ini,η_ini)

# receive solutions from all processors
ηall = MPI.Gather(η,0,comm)
if rank == 0
    ηall = reshape(ηall,Nt+1,N)
end

MPI.Barrier(comm)
MPI.Finalize()

## PLOTTING
if rank == 0
    fig,ax = subplots()

    l1, = ax[:plot](x_h,ηall[1,:])
    ax[:set_ylim](-2,2)
    ax[:set_xlabel]("x")
    ax[:set_ylabel]("η")
    tight_layout()

    for it = 2:3:Nt+1
        pause(0.00001)
        l1[:set_data](x_h,ηall[it,:])
    end

    pause(2)
    close(fig)
end
