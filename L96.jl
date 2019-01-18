function RK4(previous)
    k1 = ode(previous)
    k2 = ode(previous + 0.5*dt*k1)
    k3 = ode(previous + 0.5*dt*k2)
    k4 = ode(previous + dt*k3)

    return previous + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
end

function ode(x)
    return (circshift(x, -1) - circshift(x, 2)) .* circshift(x, 1) - x .+ F
end

const F = 8.0
const dt = 0.005

nx = 50
state = zeros(nx)
state[1] = 8.0

for i in 1:20000
    state = RK4(state)
end
