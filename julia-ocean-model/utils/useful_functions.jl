using Oceananigans
using Oceananigans.Grids: nodes, xnode, ynode
using Oceananigans.Fields: interpolate, location
using Oceananigans.Units
using Oceananigans
using Oceananigans.Grids: architecture
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils
using Printf
using CairoMakie

####
#### Interpolating the bathymetry map over a grid.
#### The bathymetry is assumed to be a 2D array of size (1440, 600)
####

function immersed_boundary_grid(grid, bathymetry)
    
    if size(bathymetry) != size(grid)
        bathymetry_grid = LatitudeLongitudeGrid(size = (1440, 600, 1),
                                                longitude = (-180, 180),
                                                latitude = (-75, 75),
                                                halo = (4, 4, 4),
                                                z = (0, 1))

        Nx, Ny, _ = size(grid)

        bathymetry_mask = zeros(Nx, Ny)

        old_field = CenterField(bathymetry_grid)
        set!(old_field, reshape(bathymetry, 1440, 600, 1))

        for i in 1:Nx, j in 1:Ny
            bathymetry_mask[i, j] = interpolate(old_field, grid.λᶜᵃᵃ[i], grid.φᵃᶜᵃ[j], grid.zᵃᵃᶜ[1])
        end
    else
        bathymetry_mask = deepcopy(bathymetry)
    end
    
    bathymetry_mask = Float64.(bathymetry_mask .>= 0)
    
    return ImmersedBoundaryGrid(grid, GridFittedBoundary(bathymetry_mask))
end

####
#### Visualization of the grid and discrete continent map
####

function visualize_cartesian_grid(grid)
    Nx, Ny, _ = size(grid)
   
    λ = zeros(Nx, Ny)
    ϕ = zeros(Nx, Ny) 
    z = ones(Nx, Ny)


    for i in 1:Ny
        λ[:, i] .= xnodes(Center, grid)
    end
    for i in 1:Nx
        ϕ[i, :] .= ynodes(Center, grid)
    end
    
    bat = Float64.(deepcopy(interior(grid.immersed_boundary.mask, :, :, 1)))
    bat[ bat .== 0 ] .= NaN
    
    fig = Figure(resolution = (1000, 500))
    
    ax1 = Axis(fig[1, 1])
    heatmap!(ax1, λ[:, 1], ϕ[1, :], bat, colormap = :hot)
    wireframe!(ax1, λ, ϕ, z, color = :black, linewidth = 0.05)
    hidedecorations!(ax1, ticklabels=false)
    return fig
end

####
#### Setting U's data to the array Um
#### Um is assumed to have the shape (1440, 600, 1) or (1440, 601, 1)
####

function set_velocity_from_array!(U, Um)
    
    grid = U.grid
    
    if any(size(U) != size(Um))
        old_grid = LatitudeLongitudeGrid(size = (1440, 600, 1),
                                         longitude = (-180, 180),
                                         latitude = (-75, 75),
                                         z = (0, 1))
        Nx, Ny, _ = size(U) 
        old_field = Field{location(U)...}(old_grid)

        tmp = zeros(size(old_field))
        
        # Allow Vm to be passed to U and Um to be passed to V
        if size(Um)[2] == size(old_field)[2] + 1
            tmp .= Um[:, 1:end-1, :]
        elseif size(Um)[2] == size(old_field)[2] - 1
            tmp[:, 1:end-1, :] .= Um
        else
            tmp .= Um
        end
        set!(old_field, tmp)
        Un = zeros(Nx, Ny, 1)

        loc = location(U)
        for i in 1:Nx, j in 1:Ny
            Un[i, j] = interpolate(old_field, xnode(i, j, 1, U), ynode(i, j, 1, U), old_grid.zᵃᵃᶜ[1])
        end
        set!(U, Un)
    else
        set!(U, Um)
    end
end

####
#### Progress function to keep track of the simulation 
####

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


#####
##### Shenanigans for time-dependent velocity simulation
#####

const vel_idx  = 362
const tot_time = 362

@inline current_time_index(time, tot_time)   = mod(unsafe_trunc(Int32, time / 10days),     tot_time) + 1
@inline next_time_index(time, tot_time)      = mod(unsafe_trunc(Int32, time / 10days) + 1, tot_time) + 1
@inline cyclic_interpolate(u₁, u₂, time)     = u₁ + mod(time / 10days, 1) * (u₂ - u₁)
    
@kernel function _set_velocities!(u_mod, v_mod, u₁, u₂, v₁, v₂, time)
    i, j, k = @index(Global, NTuple)
    u_mod[i, j, k] = cyclic_interpolate(u₁[i, j, k], u₂[i, j, k], time)
    v_mod[i, j, k] = cyclic_interpolate(v₁[i, j, k], v₂[i, j, k], time)
end

@inline function change_velocity_function!(simulation)
    time = simulation.model.clock.time

    grid = simulation.model.grid
    n₁ = current_time_index(time, tot_time)
    n₂ = next_time_index(time, tot_time)

    u₁ = ut[n₁]
    u₂ = ut[n₂]
    
    v₁ = vt[n₁]
    v₂ = vt[n₂]
    
    u_mod = simulation.model.velocities.u
    v_mod = simulation.model.velocities.v

    event = launch!(architecture(grid), grid, :xyz, _set_velocities!, u_mod, v_mod, u₁, u₂, v₁, v₂, time)
    wait(Oceananigans.Architectures.device(architecture(grid)), event)

    return nothing
end

@inline function calculate_velocities(u, v, U, V, δ)
    um = []
    vm = []

    for i in 1:length(u)
        push!(um, u[i] .* δ .+ Um .* (1 - δ))
        push!(vm, v[i] .* δ .+ Vm .* (1 - δ))
    end
    
    return um, vm
end