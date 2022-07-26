
using Oceananigans
using Oceananigans.Grids: nodes
using Oceananigans.Fields: interpolate
using Printf
using CairoMakie

function geographic2cartesian(λ, φ; r=1)
    Nλ = length(λ)
    Nφ = length(φ)

    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ) 
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

function visualize_grid(grid)

    Nx, Ny, _ = size(grid)

    λ, ϕ, z = nodes((Center, Center, Center), grid) 
    x, y, z = geographic2cartesian(λ, ϕ, r=1.01)
    
    fig = Figure(resolution = (1000, 1000))

    bat = Float64.(deepcopy(interior(grid.immersed_boundary.mask, :, :, 1)))

    bat[ bat .== 0 ] .= NaN

    range = [1:Nx, 1:Ny]
    ax1 = fig[1, 1] = LScene(fig) # make plot area wider
    wireframe!(ax1, x[range...], y[range...], z[range...], linewidth = 0.2)
    surface!(ax1, x, y, z, color = bat, colormap = :hot)
    wireframe!(ax1, λ, ϕ, z, color = :black, linewidth = 0.05)

    init = (π/5, π/3 + π/3, 0)
    rotate_cam!(ax1.scene, init)

    return fig
end

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