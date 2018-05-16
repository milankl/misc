# define constants

const g = 10.
const H = 10.
const N = 500
const Nt = 300
const cfl = 0.9
const rho = 1.
const ν = 0.006
const F0 = 500.
const Lx = 1.
const dx = Lx/N
const cph = sqrt(g*H)
const dt = cfl * dx / cph
const RKα = [1/6.,1/3.,1/3.,1/6.]
const RKβ = [0.5,0.5,1.]

# grid
const x_h = dx/2:dx:Lx
const x_u = dx:dx:Lx        # no u-point at x=0 but at x=L (periodicity)

F(x) = -F0*sin.(8π*x/Lx).*sin.(π*x/Lx)/rho
const Fx = F(x_h)

# OPERATORS
function GT(N)
    GTx = spdiagm((ones(N),-ones(N-1)),(0,-1),N,N) / dx
    GTx[1,N] = -1./dx
    GTx
end

const GTx = GT(N)
const Gux = -GTx'

const ITu = abs.(GTx*dx/2.)
const IuT = abs.(Gux*dx/2.)

# initial conditions
η = zeros(N)
u = zeros(N)

function rhs(du,dη,h_u,u,η)
    h_u = ITu*(η+H)

    # non-linear
    du = -GTx*(.5*IuT*(u.^2) + g*η - ν*Gux*u) + Fx./h_u
    dη = -Gux*(u.*h_u)

    return du,dη
end

function time_integration(Nt,u,η)

    # pre-allocate
    u0,η0 = zeros(u),zeros(η)
    u1,η1 = zeros(u),zeros(η)
    du,dη = zeros(u),zeros(η)
    h_u = zeros(u)

    for i = 1:Nt
        u1[:] = u
        η1[:] = η

        for rki = 1:4
            du,dη = rhs(du,dη,h_u,u,η)

            if rki < 4 # RHS update for the next RK-step
                u1 = u + RKβ[rki]*dt*du
                η1 = η + RKβ[rki]*dt*dη
            end

            # Summing all the RHS on the go
            u0 += RKα[rki]*dt*du
            η0 += RKα[rki]*dt*dη

        end

        u[:] = u0
        η[:] = η0
    end
    u,η
end
