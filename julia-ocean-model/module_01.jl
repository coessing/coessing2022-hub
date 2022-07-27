using Oceananigans, JLD2, CairoMakie
using Downloads: download

include("usefull_functions.jl");

# We need to define the remote path to the file
onlinepath = "https://github.com/simone-silvestri/coessing-data/raw/main/"

file_path = download(onlinepath * "bathymetry.jld2", "./bathymetry.jld2")
file_bathymetry = jldopen(file_path);

bathymetry      = file_bathymetry["bathymetry"];
bathymetry_mask = Float64.(bathymetry .>= 0);

Nlon = 1440;
Nlat = 600;

underlying_grid = LatitudeLongitudeGrid(size = (Nlon, Nlat, 1),
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (4, 4, 4),
                                        z = (0, 1));

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(bathymetry_mask));