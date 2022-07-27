using Oceananigans.Units

include("module_01.jl")

onlinepath = "https://github.com/simone-silvestri/coessing-data/raw/main/"

dh = DataDep("quarter_degree_velocity_timeseries",
     "Mean velocity fields for the global latitude longitude simulation",
     onlinepath * "prescribed_mean_fields.jld2"
)

DataDeps.register(dh)

datadep"quarter_degree_velocity_timeseries"

file_velocities = jldopen(@datadep_str "quarter_degree_velocity_timeseries/prescribed_mean_fields.jld2")

U = XFaceField(grid)
V = YFaceField(grid)

const Um = file_velocities["um"]
const Vm = file_velocities["vm"]

set!(U, Um)
set!(V, Vm)
