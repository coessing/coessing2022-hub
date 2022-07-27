# using Pkg
# Pkg.add(url="https://github.com/CliMA/Oceananigans.jl", rev="ss/particle_tracker")
# Pkg.add("DataDeps")
# Pkg.add("KernelAbstractions")
# Pkg.add("JLD2")
# Pkg.add("CairoMakie")

using Oceananigans, DataDeps, JLD2, CairoMakie


include("usefull_functions.jl")

onlinepath = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

dh = DataDep("quarter_degree_near_global_lat_lon",
     "Forcing data for global latitude longitude simulation",
     onlinepath * "bathymetry-1440x600.jld2"
)

DataDeps.register(dh)

datadep"quarter_degree_near_global_lat_lon"
file_bathymetry = jldopen(@datadep_str "quarter_degree_near_global_lat_lon/bathymetry-1440x600.jld2")

bathymetry      = file_bathymetry["bathymetry"]
bathymetry_mask = Float64.(bathymetry .>= 0)

Nx = 1440
Ny = 600

underlying_grid = LatitudeLongitudeGrid(size = (Nx, Ny, 1),
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (4, 4, 4),
                                        z = (0, 1))

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(bathymetry_mask))