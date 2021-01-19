__precompile__()

module EnhancedGJK

using StaticArrays
import GeometryTypes
const gt = GeometryTypes
using CoordinateTransformations: Transformation,
                                 transform_deriv,
                                 IdentityTransformation
using LinearAlgebra
using Statistics: mean
using Meshing
using ForwardDiff
# using DiffResults
# using GeometryBasics: Mesh, TriangleFace, GLTriangleFace, Point3f0

export CollisionCache,
       gjk!,  # signed volume method
       gjk,  # signed volume method
       body_simplex,  # return simplex points on body in body frame
       gjk_original!, # johnson method
       gjk_original,  # johnson method
       GJKResult,
       NeighborMesh,
       ReferenceDistance,
       closest_point_in_world,
       closest_point_in_body,
       separation_distance,
       simplex_penetration_distance,
       linear_combination,
       implicit_surface,
       construct_surface,
       construct_centers,
       eval_surface,
       surface_gradient!,
       fzero,
       make_mesh,
       make_geom,
       make_mesh_sdf

include("tagged_points.jl")
include("simplices.jl")
include("johnson_distance.jl")
include("neighbor_mesh.jl")
include("traits.jl")
include("gjk.jl")
include("reference_distance.jl")
include("surface.jl")

end # module
