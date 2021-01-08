using EnhancedGJK
using BenchmarkTools
using MeshIO  # v0.3.2
using FileIO
using CoordinateTransformations, StaticArrays

# Dimensions: (0.264 x 0.14 x 0.0885 m)
c1 = load(string(pwd(),"/test/meshes/r_foot_chull.obj"))
c2 = load(string(pwd(),"/test/meshes/r_foot_chull.obj"))
trans = Translation(SVector{3,Float64}(rand(3,1)))
poseA = trans
trans = Translation(SVector{3,Float64}(0,0,0))
poseB = trans

# GeometryBasics testing
result = gjk(c1, c2, poseA, poseB)
ref_result = gjk_original(c1, c2, poseA, poseB)

@btime gjk(c1, c2, poseA, poseB)
@btime gjk_original(c1, c2, poseA, poseB)

println("\n","\n")
println(result.weights, "\n", result.in_collision, "\n", result.closest_point_in_body, "\n", result.iterations, "\n", result.termination)
println(ref_result.weights, "\n", ref_result.in_collision, "\n", ref_result.closest_point_in_body, "\n", ref_result.iterations)

println(separation_distance(result))
println(separation_distance(ref_result))

println(result.mesh_simplex)
println(result.nvrtx)
