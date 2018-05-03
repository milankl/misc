
function time_integration(N,x,y,z)

    # pre-allocate
    dx,dy,dz = 0.,0.,0.

    for i = 1:N

        # RHS
        @fastmath dx = dt*(σ*(y-x))
        @fastmath dy = dt*(x*(ρ-z) - y)
        @fastmath dz = dt*(x*y - β*z)

        # Euler forward
        x += dx
        y += dy
        z += dz
    end
    x,y,z
end

N = Int(1e8)        # number of time steps
const dt = 0.01     # time step

# initial conditions
x,y,z = 5.,5.,5.

# L63 parameters
const σ = 10.
const ρ = 28.
const β = 8./3.

# compile
time_integration(1,x,y,z)

# RUN
@time x,y,z = time_integration(N,x,y,z)
println((x,y,z))
