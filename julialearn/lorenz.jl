function time_integration(N,x,y,z)

    dx,dy,dz = 0.,0.,0.         # pre-allocate

    for i = 1:N

        dx = dt*(σ*(y-x))       # RHS
        dy = dt*(x*(ρ-z) - y)
        dz = dt*(x*y - β*z)

        x += dx                 # Euler forward
        y += dy
        z += dz
    end
    x,y,z
end

N = Int(1e8)                    # number of time steps
const dt = 0.01                 # time step

x,y,z = 5.,5.,5.                # initial conditions

const σ = 10.                   # L63 parameters
const ρ = 28.
const β = 8./3.

time_integration(1,x,y,z)       # compile

@time x,y,z = time_integration(N,x,y,z) # RUN
