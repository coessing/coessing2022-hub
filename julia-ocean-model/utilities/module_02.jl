using Oceananigans.Units
using Downloads: download

include("module_01.jl")

onlinepath = "https://github.com/simone-silvestri/coessing-data/raw/main/"

file_path = download(onlinepath * "prescribed_mean_fields.jld2", "./prescribed_mean_fields.jld2")

file_velocities = jldopen(file_path)
    
U = XFaceField(grid)
V = YFaceField(grid)

const Um = file_velocities["um"]
const Vm = file_velocities["vm"]

set!(U, Um)
set!(V, Vm)
