using EnhancedGJK
using BenchmarkTools
using MeshIO  # v0.3.2
using FileIO
using CoordinateTransformations, StaticArrays
using MeshCat, Colors
import GeometryTypes: HyperRectangle, HyperSphere, Vec, Point, HomogenousMesh

# Dimensions: (0.264 x 0.14 x 0.0885 m)
c1 = load(string(pwd(),"/test/meshes/r_foot_chull.obj"))
c2 = load(string(pwd(),"/test/meshes/r_foot_chull.obj"))
# trans = Translation(SVector{3,Float64}(rand(3,1)))
trans = Translation(SVector{3,Float64}(0, 0, 0.20))
poseA = trans
trans = Translation(SVector{3,Float64}(0,0,0))
poseB = trans

# GeometryBasics testing
result = gjk(c1, c2, poseA, poseB)
ref_result = gjk_original(c1, c2, poseA, poseB)

@btime gjk(c1, c2, poseA, poseB)
@btime gjk_original(c1, c2, poseA, poseB)

@code_warntype gjk(c1, c2, poseA, poseB)
@profiler gjk(c1, c2, poseA, poseB)

println("\n","\n")
println(result.weights, "\n", result.in_collision, "\n", result.closest_point_in_body, "\n", result.iterations, "\n", result.termination)
println(ref_result.weights, "\n", ref_result.in_collision, "\n", ref_result.closest_point_in_body, "\n", ref_result.iterations)

println(separation_distance(result))
println(separation_distance(ref_result))

# Visualization
vis = Visualizer()
open(vis)

# Set objects
setobject!(vis["mesh/1"], c1)
setobject!(vis["mesh/2"], c2)
settransform!(vis["mesh/1"], poseA)
settransform!(vis["mesh/2"], poseB)

# Show distance arrows

# SV GJK
dist_vec1 = ArrowVisualizer(vis["dist/a"])
setobject!(dist_vec1, MeshLambertMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))
settransform!(dist_vec1, Point(poseA(result.closest_point_in_body.a)), Point(poseB(result.closest_point_in_body.b)))  # poseB transformation

# Reference 
dist_vec2 = ArrowVisualizer(vis["dist/b"])
setobject!(dist_vec2, MeshLambertMaterial(color=RGBA{Float32}(0, 1, 0, 0.5)))
settransform!(dist_vec2, Point(poseA(ref_result.closest_point_in_body.a)),
    Point(poseB(ref_result.closest_point_in_body.b)))  # poseB transformation
