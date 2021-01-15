using EnhancedGJK
import GeometryTypes
const gt = GeometryTypes
# import GeometryBasics
using BenchmarkTools
using MeshIO, FileIO  # v0.3.2
using CoordinateTransformations, StaticArrays
# import GeometryTypes: HyperRectangle, HyperSphere, Vec, Point, HomogenousMesh, GLNormalMesh
using MeshCat, Colors

# Dimensions: (0.264 x 0.14 x 0.0885 m)
T = Float64
c1 = load(string(pwd(),"/test/meshes/r_foot_chull.obj"))
c2 = gt.GLNormalMesh(gt.HyperSphere(gt.Point(0,0,0.1), 0.01))  # size is one-sigma sensor noise
trans = Translation(SVector{3,T}(0,0,0))
poseA = trans
trans = Translation(SVector{3,T}(0,0,0))
poseB = trans

# Calculate
result = gjk(c1, c2, poseA, poseB)
ref_result = gjk_original(c1, c2, poseA, poseB)

println(separation_distance(result))
println(separation_distance(ref_result))

# centers are constructed from the vertices of the body, since the surface is approximated
# can decide on smaller subset of vertices, or all vertices
T = Float64
env_point = SVector{1,SVector{3,T}}(SVector{3,T}(0,0,0.1))  # inertial frame
vertices = poseB.(gt.vertices(c1))  # inertial frame
signed_distances = SVector{1,T}(separation_distance(result))
centers, distances = construct_centers(vertices, env_point, signed_distances)
coeff = construct_surface(centers, distances)  # body frame

# evaluate surface and check with actual gjk results
eval_surface(coeff, centers, SVector(0,0,0.1))
eval_surface(coeff, centers, SVector(0,0,0.05))

_poseB = Translation(SVector{3,T}(0,0,-0.05))
_result = gjk(c1, c2, poseA, _poseB)
println(separation_distance(_result))

# Create implicit surface mesh
# surf_pts, surf_faces = make_geom(coeff, centers, SVector{3,T}(1,1,1), SVector{3,T}(2,2,2), 20)
# mesh = make_mesh(surf_pts, surf_faces)
mesh = make_mesh_sdf(coeff, centers, SVector{3,T}(-0.5,-0.5,-0.5), SVector{3,T}(2,2,2), 0.005)

# # get vertices, normals, and faces to create HomogenousMesh
# _vertices = GeometryBasics.decompose(GeometryBasics.Point{3,T}, mesh)
# _normals = GeometryBasics.decompose_normals(mesh)
# _faces = GeometryBasics.decompose(GeometryBasics.TriangleFace{Int}, mesh)
# _texture = [gt.Point(0.0f0,0) for i in 1:1692]
# implicit_mesh = gt.HomogenousMesh(_vertices, _normals, _faces, texture)

# Visualization
vis = Visualizer()
open(vis)

setobject!(vis["mesh/1"], c1, MeshPhongMaterial(color=RGBA{Float32}(1.0,1.0,0.9,0.9)))
setobject!(vis["mesh/2"], c2)
setobject!(vis["mesh/surface"], mesh, MeshPhongMaterial(color=RGBA{Float32}(0.0,1.0,0.0,0.9)))
# setobject!(vis["mesh/surface"], mesh)
settransform!(vis["mesh/1"], poseA)
settransform!(vis["mesh/2"], poseB)
settransform!(vis["mesh/surface"], poseB)

# Add gjk arrows
dist_vec = ArrowVisualizer(vis["dist/a"])
setobject!(dist_vec, MeshLambertMaterial(color=RGBA{Float32}(0, 1, 0, 0.5)))
settransform!(dist_vec, gt.Point(poseA(result.closest_point_in_body.a)),
    gt.Point(poseB(result.closest_point_in_body.b)))  # poseB transformation

# delete tasks
delete!(vis["mesh/surface"])
