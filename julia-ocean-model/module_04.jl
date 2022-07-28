global_path = "/Users/simonesilvestri/Dropbox/lab-nigeria/"
include(global_path * "module_03.jl")

using Statistics
using JLD2
using Printf
using CairoMakie
using Oceananigans.Units   
using Oceananigans.Architectures: arch_array
using Oceananigans.Operators

#####
##### Include Particles
#####

using Oceananigans.Models.HydrostaticFreeSurfaceModels: OnlyParticleTrackingModel

function run_particle_simulation!(λ₀, φ₀, n_particles, degree_spread_λ, degree_spread_φ, δ_turb, Nyears)
    
    ##### 
    ##### Setting particle's initial position
    #####

    λₚ = λ₀ .+ degree_spread_λ .* (rand(n_particles) .- 0.5);
    φₚ = φ₀ .+ degree_spread_φ .* (rand(n_particles) .- 0.5);
 
    zₚ = 0.5 .* ones(n_particles);
    lagrangian_particles = LagrangianParticles(x=λₚ, y=φₚ, z=zₚ)

    #####
    ##### Uploading velocity depending on δ_turb (share of mean vs instantaneous velocity)
    #####

    if δ_turb !=1 && δ_turb !=0
        for i in 1:vel_idx
            @info "$i of $vel_idx"
            un[i] = δ_turb .* u[i] .+ (1 - δ_turb) .* Um
            vn[i] = δ_turb .* v[i] .+ (1 - δ_turb) .* Vm
        end
    end

    #####
    ##### Setup model
    #####

    model = HydrostaticFreeSurfaceModel(grid = grid,
                                        velocities = PrescribedVelocityFields(u = U, v = V),
                                        coriolis = nothing,
                                        buoyancy = nothing,
                                        closure = nothing,
                                        tracers = (),
                                        particles = lagrangian_particles)
    @info "model initialized"
    @info "model only advects particles $(model isa OnlyParticleTrackingModel)"

    #####
    ##### Simulation setup
    #####

    Δt = 6hours

    simulation = Simulation(model, Δt = Δt, stop_time = Nyears*years)

    start_time = [time_ns()]

    # If we use a velocity which is a function of time, remember to change it at every timestep!
    if δ_turb != 0 
        simulation.callbacks[:interp_vel] = Callback(velocity_function, IterationInterval(1))
    end

    # Callback to keep track of the progress
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

    save_interval = 10days

    simulation.output_writers[:surface_fields] = JLD2OutputWriter(model,  (; u = model.velocities.u, v = model.velocities.v, particles=model.particles,),
                                                                schedule = TimeInterval(save_interval),
                                                                filename = "output_particles",
                                                                overwrite_existing = true)

    # Let's goo!
    @info "Running with Δt = $(prettytime(simulation.Δt))"

    run!(simulation)

    @info """

        Simulation took $(prettytime(simulation.run_wall_time))
        Time step: $(prettytime(Δt))
    """
    return model
end

####
#### Plotting results
####

function visualize_results(output_prefix)
    bat = deepcopy(arch_array(CPU(), bathymetry_mask))

    surface_file = jldopen(output_prefix * ".jld2")

    bat = Float64.(bat)
    bat[bat .== 1] .= NaN

    iterations = parse.(Int, keys(surface_file["timeseries/t"]))

    iter = Observable(0)

    xi(iter) = surface_file["timeseries/particles/" * string(iter)].x
    yi(iter) = surface_file["timeseries/particles/" * string(iter)].y
    ti(iter) = string(surface_file["timeseries/t/" * string(iter)] / day)

    ui(iter) = surface_file["timeseries/u/" * string(iter)][:, :,       1]
    vi(iter) = surface_file["timeseries/v/" * string(iter)][:, 1:end-1, 1]
    
    Px = @lift Array(xi($iter))
    Py = @lift Array(yi($iter))

    x₀ = []
    y₀ = []
    for j in eachindex(iterations)
        push!(x₀, xi(iterations[j]))
        push!(y₀, yi(iterations[j]))
    end

    speed = @lift (Array(ui($iter)).^2 .+ Array(vi($iter)).^2).^(0.5) .+ bat

    fig = Figure(resolution = (1400, 900))

    λ = range(-179.75, 179.75, length = 1440)
    φ = range(-74.75, 74.75, length = 600)

    ax = Axis(fig[1, 1], title="Tracer concentration (m)")
    hm = CairoMakie.heatmap!(ax, λ, φ, speed, colorrange=(0.0, 0.5), colormap = :viridis, nan_color = :black, interpolate = true)
    CairoMakie.scatter!(ax, Px, Py, color = :red)

    CairoMakie.record(fig, output_prefix * ".mp4", 1:length(iterations), framerate=8) do i
        @info "Plotting iteration $i of $(length(iterations))..."
        if i > 1
            CairoMakie.scatter!(ax, x₀[i-1], y₀[i-1], color = :white, markersize = 4)
        end
        iter[] = iterations[i]
    end

    display(fig)

    close(surface_file)
end
