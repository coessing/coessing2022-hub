using Oceananigans
using Oceananigans.Grids: nodes
using Oceananigans.Fields: interpolate
using Printf
using CairoMakie

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
    hidedecorations!(ax1)

    return fig
end

function immersed_boundary_grid(grid, bathymetry)
    
    bat_grid = LatitudeLongitudeGrid(size = (1440, 600, 1),
                                     longitude = (-180, 180),
                                     latitude = (-75, 75),
                                     halo = (4, 4, 4),
                                     z = (0, 1))

    Nx, Ny, _ = size(grid)

    bat_new = zeros(Nx, Ny)

    old_field = Field{Center, Center, Center}(bat_grid)
    set!(old_field, reshape(bathymetry, 1440, 600, 1))
    
    for i in 1:Nx, j in 1:Ny
        bat_new[i, j] = interpolate(old_field, grid.λᶜᵃᵃ[i], grid.φᵃᶜᵃ[j], bat_grid.zᵃᵃᶜ[1])
    end

    bat_new = Float64.(bat_new .>= 0)

    return ImmersedBoundaryGrid(grid, GridFittedBoundary(bat_new))
end
