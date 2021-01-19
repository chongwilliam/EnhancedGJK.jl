using EnhancedGJK
import GeometryTypes
const gt = GeometryTypes
using BenchmarkTools
using MeshIO, FileIO  # v0.3.2
using CoordinateTransformations, StaticArrays
using MeshCat, Colors

# Dimensions: (0.264 x 0.14 x 0.0885 m)
T = Float64
c = []
x, y, z = 0.05, 0.05, 0.01
push!(c, load(string(pwd(),"/test/meshes/r_foot_chull.obj")))
push!(c, gt.GLNormalMesh(gt.HyperSphere(gt.Point(0.,0,z), 0.005)))  # size is one-sigma sensor noise
poseA = Translation(SVector{3,T}(0,0,0))
poseB = Translation(SVector{3,T}(0,0,0))

# Calculate
result = []
push!(result, gjk(c[1], c[2], poseA, poseB))

# Additional points
push!(c, gt.GLNormalMesh(gt.HyperSphere(gt.Point(0,y,z), 0.005)))
push!(result, gjk(c[1], c[3], poseA, poseB))

push!(c, gt.GLNormalMesh(gt.HyperSphere(gt.Point(x,y,z), 0.005)))
push!(result, gjk(c[1], c[4], poseA, poseB))

# Construct surfaces
T = Float64
env_point = SVector{3,SVector{3,T}}(SVector{3,T}(0,0,z), SVector{3,T}(0,y,z), SVector{3,T}(x,y,z))  # inertial frame
vertices = poseA.(gt.vertices(c[1]))  # inertial frame
signed_distances = SVector{3,T}(result[1].signed_distance,
                                result[2].signed_distance,
                                result[3].signed_distance)
coeff, centers = implicit_surface(vertices, env_point, signed_distances)
# centers, distances = construct_centers(vertices, env_point, signed_distances)
# coeff = construct_surface(centers, distances)  # body frame

# Benchmark
@btime construct_centers(vertices, env_point, signed_distances)
@btime construct_surface(centers, distances)
@btime implicit_surface(vertices, env_point, signed_distances)

# Evaluate surface and check with actual gjk results
eval_surface(coeff, centers, env_point[1])

# Differentiation
pt = [0, 0, 0.01]
∇f = similar(pt)
surface_gradient!(∇f, coeff, centers, pt)
print(∇f)
# ∇f, H = surface_gradient(coeff, centers, SVector{3,T}(0,0,0.1))
# @btime surface_gradient(coeff, centers, SVector{3,T}(0,0,0.1))

# Test intersection set algorithm with sampled points around points in collision
zero_set, cnt = fzero(coeff, centers, [env_point[1]+[0,0,0.01]], ∇f)
cnt
eval_surface(coeff, centers, zero_set[1])
# eval_surface(coeff, centers, zero_set[2])
# eval_surface(coeff, centers, zero_set[3])


# Create implicit surface mesh
# surf_pts, surf_faces = make_geom(coeff, centers, SVector{3,T}(1,1,1), SVector{3,T}(2,2,2), 20)
# mesh = make_mesh(surf_pts, surf_faces)
mesh = make_mesh_sdf(coeff, centers, SVector{3,T}(-0.5,-0.5,-0.5), SVector{3,T}(2,2,2), 0.01)

# # get vertices, normals, and faces to create HomogenousMesh
# _vertices = GeometryBasics.decompose(GeometryBasics.Point{3,T}, mesh)
# _normals = GeometryBasics.decompose_normals(mesh)
# _faces = GeometryBasics.decompose(GeometryBasics.TriangleFace{Int}, mesh)
# _texture = [gt.Point(0.0f0,0) for i in 1:1692]
# implicit_mesh = gt.HomogenousMesh(_vertices, _normals, _faces, texture)

# Visualization
vis = Visualizer()
open(vis)

setobject!(vis["mesh/1"], c[1], MeshPhongMaterial(color=RGBA{Float32}(1.0,1.0,0.9,0.9)))
setobject!(vis["mesh/2"], c[2])
setobject!(vis["mesh/3"], c[3])
setobject!(vis["mesh/4"], c[4])
setobject!(vis["mesh/surface"], mesh, MeshPhongMaterial(color=RGBA{Float32}(0.0,0.0,1.0,0.9)))

settransform!(vis["mesh/1"], poseA)
settransform!(vis["mesh/2"], poseB)
settransform!(vis["mesh/3"], poseB)
settransform!(vis["mesh/4"], poseB)
settransform!(vis["mesh/surface"], poseA)

# Add gjk arrows
dist_vec = ArrowVisualizer(vis["dist/a"])
setobject!(dist_vec, MeshLambertMaterial(color=RGBA{Float32}(0,1,0,0.5)))
settransform!(dist_vec, gt.Point(poseA(result.closest_point_in_body.a)),
    gt.Point(poseB(result.closest_point_in_body.b)))  # poseB transformation

# delete tasks
delete!(vis["mesh/1"])
delete!(vis["mesh/2"])
delete!(vis["mesh/3"])
delete!(vis["mesh/4"])
delete!(vis["mesh/surface"])

# delete visualizer
delete!(vis)
