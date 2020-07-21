using LinearAlgebra
using Convex
using ECOS
using BSON: @load

@load "test/msl_traj.bson" x_traj u_traj t_traj
