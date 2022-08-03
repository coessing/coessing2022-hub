function visualize_mean_velocity_results(model)
  
    grid = model.grid 
    bat  = grid.immersed_boundary.mask[3:end-1, 3:end-1, 3]

    bat = Float64.(bat)
    bat[bat .== 1] .= NaN
    
    iterations = length(particles_output)
        
    xi(iter) = particles_output[iter].properties.x
    yi(iter) = particles_output[iter].properties.y

    Nx, Ny, Nz = size(grid)
    
    u = interior(model.velocities.u, :, 1:Ny, 1)
    v = interior(model.velocities.v, :, 1:Ny, 1)
    
    iter = Observable(1)
    
    Px = @lift Array(xi($iter))
    Py = @lift Array(yi($iter))

    x₀ = []
    y₀ = []
    for j in 1:length(iterations)
        push!(x₀, xi(j))
        push!(y₀, yi(j))
    end
    
    Nt = 20

    Px_vec = Vector(undef, Nt)
    Py_vec = Vector(undef, Nt)

    for j in 1:Nt
        Px_vec[j] = @lift $iter - j * 3 > 1 ? Array(xi($iter - j * 3)) : Array(xi(1))
        Py_vec[j] = @lift $iter - j * 3 > 1 ? Array(yi($iter - j * 3)) : Array(yi(1))
    end
    
    speed = 0.5 .* (u.^2 .+ v.^2).^(0.5) .+ bat

    fig = Figure(resolution = (1000, 500))
    
    λ = range(-179.75, 179.75, length = Nx)
    φ = range(-74.75, 74.75, length = Ny)

    ax = Axis(fig[1, 1], title="Particle paths - Contour speed (√u² + v²) ms⁻¹")
    hm = heatmap!(ax, λ, φ, speed, colorrange=(0.0, 0.2), colormap = :viridis, nan_color = :black, interpolate = true)
    scatter!(ax, Px, Py, color = :red, markersize = 7)
    for j in 1:Nt
        scatter!(ax, Px_vec[j], Py_vec[j], color = (:white, (Nt - j) / Nt), markersize = 4 * (Nt - j) / Nt)
    end

    record(fig, "output_particles.mp4", 1:iterations, framerate=12) do i
        if mod(i, 50) == 0
            @info "Plotting iteration $i of $(iterations)..."
        end
        iter[] = i
    end

    @info "video recorded in output_particles.mp4"
end


function visualize_results(grid)
    
    bat = grid.immersed_boundary.mask[3:end-1, 3:end-1, 3]

    bat = Float64.(bat)
    bat[bat .== 1] .= NaN
    
    iterations = length(u_output)
        
    xi(iter) = particles_output[iter].properties.x
    yi(iter) = particles_output[iter].properties.y

    Nx, Ny, Nz = size(u_output[1])
    
    ui(iter) = interior(u_output[iter], :, 1:Ny, 1)
    vi(iter) = interior(v_output[iter], :, 1:Ny, 1)
    
    iter = Observable(1)
    
    Px = @lift Array(xi($iter))
    Py = @lift Array(yi($iter))

    x₀ = []
    y₀ = []
    for j in 1:length(iterations)
        push!(x₀, xi(j))
        push!(y₀, yi(j))
    end
    
    Nt = 20

    Px_vec = Vector(undef, Nt)
    Py_vec = Vector(undef, Nt)

    for j in 1:Nt
        Px_vec[j] = @lift $iter - j * 3 > 1 ? Array(xi($iter - j * 3)) : Array(xi(1))
        Py_vec[j] = @lift $iter - j * 3 > 1 ? Array(yi($iter - j * 3)) : Array(yi(1))
    end
    
    speed = @lift (Array(ui($iter)).^2 .+ Array(vi($iter)).^2).^(0.5) .+ bat

    fig = Figure(resolution = (1000, 500))
    
    λ = range(-179.75, 179.75, length = Nx)
    φ = range(-74.75, 74.75, length = Ny)

    ax = Axis(fig[1, 1], title="Particle paths - Contour speed (√u² + v²) ms⁻¹")
    hm = heatmap!(ax, λ, φ, speed, colorrange=(0.0, 0.5), colormap = :viridis, nan_color = :black, interpolate = true)
    scatter!(ax, Px, Py, color = :red, markersize = 7)
    for j in 1:Nt
        scatter!(ax, Px_vec[j], Py_vec[j], color = (:white, (Nt - j) / Nt), markersize = 4 * (Nt - j) / Nt)
    end

    record(fig, "output_particles.mp4", 1:iterations, framerate=12) do i
        if mod(i, 50) == 0
            @info "Plotting iteration $i of $(iterations)..."
        end
        iter[] = i
    end

    @info "video recorded in output_particles.mp4"
end

using Base64: base64encode

function display_mp4(filename)
    display("text/html", string("""<video autoplay controls><source src="data:video/x-m4v;base64,""",
    base64encode(open(read,filename)),"""" type="video/mp4"></video>"""))
end
