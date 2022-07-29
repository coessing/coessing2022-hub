
function visualize_results(output_prefix)
    
    surface_file = jldopen(output_prefix * ".jld2")

    bat = surface_file["grid/immersed_boundary/mask"][4:end-3, 4:end-3, 4]

    bat = Float64.(bat)
    bat[bat .== 1] .= NaN

    iterations = parse.(Int, keys(surface_file["timeseries/t"]))

    iter = Observable(1)

    xi(iter) = surface_file["timeseries/particles/" * string(iterations[iter])].x
    yi(iter) = surface_file["timeseries/particles/" * string(iterations[iter])].y
    ti(iter) = string(surface_file["timeseries/t/" * string(iterations[iter])] / day)

    ui(iter) = surface_file["timeseries/u/" * string(iterations[iter])][:, :,       1]
    vi(iter) = surface_file["timeseries/v/" * string(iterations[iter])][:, 1:end-1, 1]
    
    Px = @lift Array(xi($iter))
    Py = @lift Array(yi($iter))

    x₀ = []
    y₀ = []
    for j in 1:length(iterations)
        push!(x₀, xi(j))
        push!(y₀, yi(j))
    end

    Px_vec = Vector(undef, 10)
    Py_vec = Vector(undef, 10)
    for j in 1:10
        Px_vec[j] = @lift $iter - j * 3 > 1 ? Array(xi($iter - j * 3)) : Array(xi(1))
        Py_vec[j] = @lift $iter - j * 3 > 1 ? Array(yi($iter - j * 3)) : Array(yi(1))
    end


    speed = @lift (Array(ui($iter)).^2 .+ Array(vi($iter)).^2).^(0.5) .+ bat

    fig = Figure(resolution = (1000, 500))

    λ = range(-179.75, 179.75, length = 1440)
    φ = range(-74.75, 74.75, length = 600)

    ax = Axis(fig[1, 1], title="Particle paths - Contour speed (√u² + v²) ms⁻¹")
    hm = heatmap!(ax, λ, φ, speed, colorrange=(0.0, 0.5), colormap = :viridis, nan_color = :black, interpolate = true)
    scatter!(ax, Px, Py, color = :yellow, markersize = 4)
    for j in 1:10
        scatter!(ax, Px_vec[j], Py_vec[j], color = :white, markersize = 2)
    end    
    
    record(fig, output_prefix * ".mp4", 1:length(iterations), framerate=16) do i
        if mod(i, 50) == 0 
            @info "Plotting iteration $i of $(length(iterations))..."
        end
        iter[] = i
    end

    close(surface_file)
    
    @info "filename = $(output_prefix).mp4"
end

using Base64: base64encode

function display_mp4(filename)
    display("text/html", string("""<video autoplay controls><source src="data:video/x-m4v;base64,""",
    base64encode(open(read,filename)),"""" type="video/mp4"></video>"""))
end
