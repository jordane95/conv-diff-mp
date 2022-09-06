using Plots

nx = 81 # grid numbers
ny = 81
dx = 1 / (nx - 1) # real length = 2
dy = 1 / (ny - 1)
nt = 300 # total time steps
dt = .001
c = 1 # velocity
a = -3
b = 3
k = 0.01

# parameters for numerical solution
k1 = k*dt/dx^2
k2 = k*dt/dy^2
k3 = a*dt/dx
k4 = b*dt/dy

# for visualization
crange = (0, 3)

# initial condition
u = ones(nx, ny)
u[floor(Int, .4/dx):floor(Int, .6/dx + 1), floor(Int, .4/dy):floor(Int, .6/dy + 1)] .= 2 # set 0.5-1 to 2

heatmap(u, clim=crange)

@gif for n in 1:nt
    un = copy(u) # at time step n
    for i in 1:nx, j in 1:ny
        # periodic boundary condition
        lx = i-1; rx = i+1
        ly = j-1; ry = j+1
        if lx < 1; lx += nx; end
        if ly < 1; ly += ny; end
        if rx > nx; rx -= nx; end
        if ry > ny; ry -= ny; end
        # adapte for negative a, b
        u[i, j] = k1*(un[rx, j]+un[lx, j]) + k2*(un[i, ry]+un[i, ly]) + (1-2*k1-2*k2-abs(k3)-abs(k4))*un[i, j]
        if k3 >= 0
            u[i, j] += abs(k3)*un[lx, j]
        else
            u[i, j] += abs(k3)*un[rx, j]
        end
        if k4 >= 0
            u[i, j] += abs(k4)*un[i, ly]
        else
            u[i, j] += abs(k4)*un[i, ry]
        end
    end
    heatmap(u, clim=crange) # at time step n+1
end
