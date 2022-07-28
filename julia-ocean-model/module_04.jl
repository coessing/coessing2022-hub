include("module_03.jl")

using Statistics
using JLD2
using Printf
using CairoMakie
using Oceananigans.Units   
using Oceananigans.Architectures: arch_array
using Oceananigans.Operators

## Load in variable velocity fields (the have to be locally downloaded from 
## https://drive.google.com/file/d/196laviLFLe1YDQetj4oZRVO9jAMCZwjv/view?usp=sharing
## 2.3 GB of data!!!!!

file_turb_vel = jldopen("prescribed_surface_velocities.jl")

const u = file_turb_vel["u"]
const v = file_turb_vel["v"]

const vel_idx  = length(u)
const tot_time = length(u) 

un = deepcopy(u)
vn = deepcopy(v)

#####
##### Include Particles
#####

include("visualize_particles.jl")

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
    return nothing
end

λ₀ = 55.2
φ₀ = 8.3

degree_spread_λ = 5.0
degree_spread_φ = 5.0

n_particles = 50

δ_turb = 0.8

run_particle_simulation!(λ₀, φ₀, n_particles, degree_spread_λ, degree_spread_φ, δ_turb, 10)

visualize_results("output_particles")