using Oceananigans, DataDeps, JLD2, CairoMakie

include("usefull_functions.jl");

onlinepath = "https://github.com/simone-silvestri/coessing-data/raw/";

dh = DataDep("quarter_degree_bathymetry",
     "Bathymetry for global latitude longitude simulation",
     onlinepath * "bathymetry.jld2"
);

DataDeps.register(dh);

datadep"quarter_degree_bathymetry";
file_bathymetry = jldopen(@datadep_str "quarter_degree_bathymetry/bathymetry.jld2");

bathymetry      = file_bathymetry["bathymetry"];
bathymetry_mask = Float64.(bathymetry .>= 0);

Nx = 1440;
Ny = 600;

underlying_grid = LatitudeLongitudeGrid(size = (Nx, Ny, 1),
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (4, 4, 4),
                                        z = (0, 1));

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(bathymetry_mask));