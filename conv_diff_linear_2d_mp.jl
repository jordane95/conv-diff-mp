using Distributed
addprocs(4)

@everywhere using SharedArrays

using Plots

@everywhere nx = 81
@everywhere ny = 81
@everywhere dx = 1 / (nx - 1) # space
@everywhere dy = 1 / (ny - 1)
@everywhere nt = 100
@everywhere dt = .001 # time
@everywhere a = 1
@everywhere b = 1 # velocity
@everywhere k = 0.01 # visocity

# coefficients for numerical recurrence
@everywhere k1 = k*dt/dx^2
@everywhere k2 = k*dt/dy^2
@everywhere k3 = a*dt/dx
@everywhere k4 = b*dt/dy


# get a range of array for each proc
@everywhere function myrange(q::SharedArray)
    idx = indexpids(q) # current process id
    if idx == 0
        return 1:0
    end
    nchunks = length(procs(q)) # number of procs
    splits = [round(Int, s) for s in range(0, stop=size(q,1), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end


@everywhere function advection_chunk!(un, u, irange, jrange)
    # @show (irange)  # display so we can see what's happening
    for i in irange, j in jrange
        # if i-1>=1 && i+1<=nx && j-1>=1 && j+1<=ny
        #     u[i, j] = k1*un[i+1, j] + k2*un[i, j+1] + (k1+k3)*un[i-1, j] + (k2+k4)*un[i, j-1] + (1-2*k1-2*k2-k3-k4)*un[i, j]
        # end

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
end


@everywhere advection_shared_chunk!(un, u) = advection_chunk!(un, u, myrange(u), 1:size(u,2))


# sync in time dimension, async in space dimension
function advection_shared!(un, u)
    @sync for n in 1:nt
        for p in procs(un)
            @async remotecall_wait(advection_shared_chunk!, p, un, u)
        end
        un .= u
    end
    return un
end

function advection_serial!(un, u)
    for n in 1:nt
        advection_chunk!(un, u, 1:size(u,1), 1:size(u,2));
        un .= u
    end
    return un
end

# set x=0.4-0.6, y=0.4-0.6 to 2 otherwise 1
init = u -> (u .= 1.0; u[floor(Int, .4/dx):floor(Int, .6/dx + 1), floor(Int, .4/dy):floor(Int, .6/dy + 1)] .= 2.0)

# two shared array to store current and next time step fluid filed
un = SharedArray{Float64, 2}((nx, ny), init=init) # old array
u = SharedArray{Float64, 2}((nx, ny), init=init) # new array


@time u_mp = advection_shared!(un,u);

un = SharedArray{Float64, 2}((nx, ny), init=init) # old array
u = SharedArray{Float64, 2}((nx, ny), init=init) # new array

@time u_sp = advection_serial!(un,u);

eps = 0.0
for (mp, sp) in zip(u_mp, u_sp)
    global eps += (mp-sp)^2
end
println(eps / (nx*ny))
