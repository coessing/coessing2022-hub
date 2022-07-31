using Oceananigans
using Oceananigans.Units
using Oceananigans.Diagnostics: accurate_cell_advection_timescale

include("utils/useful_functions.jl")

using JLD2
file_bathymetry = jldopen("bathymetry.jld2")

bathymetry = file_bathymetry["bathymetry"]
Nlon = 1440
Nlat = 600

initial_grid = LatitudeLongitudeGrid(size = (Nlon, Nlat, 1), longitude = (-180, 180), latitude = (-75, 75), z = (0, 1))
grid = immersed_boundary_grid(initial_grid, bathymetry);

file_velocities = jldopen("prescribed_mean_fields.jld2")

const Um = file_velocities["um"]
const Vm = file_velocities["vm"]

# By convention, U velocity is positioned at x `Face`s and `V` at y `Face`s
u_location = (Face, Center, Center)
v_location = (Center, Face, Center)
w_location = (Center, Center, Face)

U = Field(u_location, grid)
V = Field(v_location, grid)
W = Field(w_location, grid)

set_velocity_from_array!(U, Um)
set_velocity_from_array!(V, Vm)