using Oceananigans.Units
global_path = "/Users/simonesilvestri/Dropbox/lab-nigeria/"

include(global_path * "module_01.jl")

file_velocities = jldopen(global_path * "prescribed_surface_fields.jld2")

U = XFaceField(grid)
V = YFaceField(grid)

const Um = file_velocities["um"]
const Vm = file_velocities["vm"]

const u = file_velocities["u"]
const v = file_velocities["v"]

const vel_idx  = length(u)
const tot_time = length(u) 

set!(U, Um)
set!(V, Vm)

un = deepcopy(u)
vn = deepcopy(v)