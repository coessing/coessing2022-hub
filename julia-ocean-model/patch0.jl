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

function run_particle_simulation!(λ₀, φ₀, n_particles, degree_spread_λ, degree_spread_φ, δ_turb, Nyears, outputs)
    
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
    ##### Time step the simulation
    #####

    Δt = 6hours

    println("$outputs")

    simulation = Simulation(model; Δt, stop_time = Nyears * years)

    if δ_turb != 0 
        simulation.callbacks[:update_velocities] = Callback(velocity_function!, IterationInterval(1))
    end

    simulation.callbacks[:println] = Callback(write_output!, TimeInterval(outputs*days))

    run!(simulation)

    println("-1, $n_particles")
end

@inline function write_output!(simulation)
    model = simulation.model
    n_particles = length(model.particles)

    time_in_days = model.clock.time / 3600 / 24
    println("$time_in_days, $n_particles")

    positions = model.particles.properties

    for p in 1:n_particles
        println("$(positions.x[p]), $(positions.y[p])")
    end
end

# Read cmd line arguments

depth = ARGS[1]
n_particles            = parse(Int, ARGS[2])
spread, φ₀, λ₀, δ_turb = parse.(Float64, ARGS[3:6])

run_particle_simulation!(λ₀, φ₀, n_particles, spread, spread, δ_turb, 10, 5)