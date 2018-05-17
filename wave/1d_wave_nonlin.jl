# define constants
using PyPlot

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
const dxinv = 1./dx
const cph = sqrt(g*H)
const dt = cfl * dx / cph
const RKα = [1/6.,1/3.,1/3.,1/6.]
const RKβ = [0.5,0.5,1.]

# grid
const x_h = dx/2:dx:Lx
const x_u = dx:dx:Lx        # no u-point at x=0 but at x=L (periodicity)

F(x,t) = -F0*sin.(8π*x/Lx + 200*t).*sin.(π*x/Lx).^2/rho

# OPERATORS
function Gux(res,u)
    res[1] = dxinv*(u[1]-u[end])
    res[2:end] = dxinv*(u[2:end]-u[1:end-1])
    return res
end

function GTx(res,η)
    res[1:end-1] = dxinv*(η[2:end]-η[1:end-1])
    res[end] = dxinv*(η[1]-η[end])
    return res
end

function ITu(res,η)
    res[1:end-1] = .5*(η[1:end-1] + η[2:end])
    res[end] = .5*(η[1]+η[end])
    return res
end

function IuT(res,u)
    res[1] = .5*(u[1] + u[end])
    res[2:end] = .5*(u[1:end-1]+u[2:end])
    return res
end


# initial conditions
η0 = zeros(N)
u0 = zeros(N)

function rhs(du,dη,h_u,u_h,dudx,u,η,t)
    h_u = ITu(h_u,η+H)
    u_h = IuT(u_h,u.^2)
    dudx = Gux(dudx,u)

    du = -GTx(du,.5*u_h + g*η - ν*dudx)
    du += F(x_h,t)./h_u

    dη = -Gux(dη,u.*h_u)

    return du,dη
end

function time_integration(Nt,u,η)

    # pre-allocate memory
    u0,η0 = zeros(u),zeros(η)
    u1,η1 = zeros(u),zeros(η)
    du,dη = zeros(u),zeros(η)
    h_u = zeros(u)
    u_h = zeros(η)
    dudx = zeros(η)
    t = 0.

    # for output
    u_out = zeros(Nt+1,length(u))
    η_out = zeros(Nt+1,length(η))

    # store the initial conditions
    u_out[1,:] = u
    η_out[1,:] = η

    for i = 1:Nt
        u1[:] = u
        η1[:] = η

        for RKi = 1:4
            du,dη = rhs(du,dη,h_u,u_h,dudx,u1,η1,t)

            if RKi < 4 # RHS update for the next RK-step
                u1 = u + RKβ[RKi]*dt*du
                η1 = η + RKβ[RKi]*dt*dη
            end

            # Summing all the RHS on the go
            u0 += RKα[RKi]*dt*du
            η0 += RKα[RKi]*dt*dη

        end

        u[:] = u0
        η[:] = η0
        t += dt

        # # output
        u_out[i+1,:] = u
        η_out[i+1,:] = η

    end
    u_out,η_out
end

time_integration(1,u0,η0)
@time u,η = time_integration(Nt,u0,η0)

## PLOTTING
fig,ax = subplots()

l1, = ax[:plot](x_h,η[1,:])
ax[:set_ylim](-2,2)
ax[:set_xlabel]("x")
ax[:set_ylabel]("η")
tight_layout()

for it = 2:5:Nt+1
    pause(0.00001)
    l1[:set_data](x_h,η[it,:])
end

pause(2)
close(fig)
