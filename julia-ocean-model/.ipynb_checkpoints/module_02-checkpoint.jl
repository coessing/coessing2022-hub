using Oceananigans.Units

include("module_01.jl")

# Figure out where to store this data?
onlinepath = "https://github.com/simone-silvestri/coessing-data"

dh = DataDep("quarter_degree_velocity_timeseries",
     "global latitude longitude simulation: 10 year velocity timeseries",
     onlinepath * "bathymetry-1440x600.jld2"
)

file_velocities = jldopen("prescribed_surface_fields.jld2")

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