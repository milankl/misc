function rhs(dx,dy,dz,x,y,z)
    dx = σ*(y-x)
    dy = x*(ρ-z) - y
    dz = x*y - β*z
    return dx,dy,dz
end

function time_integration(N,x,y,z)
    dx,dy,dz = zeros((x,y,z))

    for i = 1:N
        dx,dy,dz = rhs(dx,dy,dz,x,y,z)
        x += dt*dx
        y += dt*dy
        z += dt*dz
    end
    return x,y,z
end

N = Int(1e8)
x,y,z = 5.,5.,5.

const σ = 10.
const ρ = 28.
const β = 8./3.

const dt = 0.01

time_integration(1,x,y,z)

@time x,y,z = time_integration(N,x,y,z)
println((x,y,z))
