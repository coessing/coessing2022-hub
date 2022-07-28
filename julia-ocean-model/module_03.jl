using Oceananigans
using Oceananigans.Grids: architecture
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils
global_path = "/Users/simonesilvestri/Dropbox/lab-nigeria/"
include(global_path * "module_02.jl")

start_time = [time_ns()]

@inline function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    @info @sprintf("Time: % 12s, iteration: %d, wall time: %s",
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration,
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

@inline current_time_index(time, tot_time)   = mod(unsafe_trunc(Int32, time / 10days),     tot_time) + 1
@inline next_time_index(time, tot_time)      = mod(unsafe_trunc(Int32, time / 10days) + 1, tot_time) + 1
@inline cyclic_interpolate(u₁, u₂, time)     = u₁ + mod(time / 10days, 1) * (u₂ - u₁)
    
@kernel function _set_velocities!(u_mod, v_mod, u₁, u₂, v₁, v₂, time)
    i, j, k = @index(Global, NTuple)
    u_mod[i, j, k] = cyclic_interpolate(u₁[i, j, k], u₂[i, j, k], time)
    v_mod[i, j, k] = cyclic_interpolate(v₁[i, j, k], v₂[i, j, k], time)
end

@inline function velocity_function!(simulation)
    time = simulation.model.clock.time

    grid = simulation.model.grid
    n₁ = current_time_index(time, tot_time)
    n₂ = next_time_index(time, tot_time)

    u₁ = un[n₁]
    u₂ = un[n₂]
    
    v₁ = vn[n₁]
    v₂ = vn[n₂]
    
    u_mod = simulation.model.velocities.u
    v_mod = simulation.model.velocities.v

    event = launch!(architecture(grid), grid, :xyz, _set_velocities!, u_mod, v_mod, u₁, u₂, v₁, v₂, time)
    wait(Oceananigans.Architectures.device(architecture(grid)), event)

    return nothing
end
