using Plots

nx = 81 # grid numbers
dx = 2 / (nx - 1) # real length = 2
nt = 25 # total time steps
dt = .025
c = 1 # velocity

u = ones(nx) # set all to one at t=0

u[floor(Int, .5/dx):floor(Int, 1/dx + 1)] .= 2 # set 0.5-1 to 2

print(u)
plot(u)

@gif for n in 1:nt
    un = copy(u) # at time step n
    for i in 2:nx
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
    end
    plot(u) # at time step n+1
end
