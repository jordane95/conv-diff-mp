using Distributed
addprocs(4)

@everywhere using SharedArrays

using Plots

@everywhere nx = 81 # grid numbers
@everywhere dx = 2 / (nx - 1) # real length = 2
@everywhere nt = 25 # total time step
@everywhere dt = .025
@everywhere c = 1 # velocity

init = u -> (
    u .= 1.0;
    u[floor(Int, .5/dx):floor(Int, 1/dx + 1)] .= 2.0
) # set 0.5-1 to 2 otherwise to 1

# two shared array to store current and next time step fluid filed
un = SharedArray{Float64}(nx, init=init) # old array
u = SharedArray{Float64}(nx, init=init) # new array


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


@everywhere function advection_chunk!(un, u, irange)
    # @show (irange)  # display so we can see what's happening
    for i in irange
        if i-1>=1 && i<=nx
            u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
        end
    end
end


@everywhere advection_shared_chunk!(un, u) = advection_chunk!(un, u, myrange(u))


# sync in time dimension, async in space dimension
function advection_shared!(un, u)
    @sync for n in 1:nt
        for p in procs(un)
            @async remotecall_wait(advection_shared_chunk!, p, un, u)
        end
        un .= u
    end
end;

function advection_serial!(un, u)
    advection_chunk!(un, u, 1:size(u,1));
end;


# @time advection_shared!(un,u);

@time advection_serial!(un,u);