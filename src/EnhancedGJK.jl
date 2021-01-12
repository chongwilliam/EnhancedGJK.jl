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

export CollisionCache,
       gjk!,  # signed volume method
       gjk,  # signed volume method
       body_simplex,  # return simplex points on body
       gjk_original,  # johnson method
       gjk_original!, # johnson method
       GJKResult,
       GJKParams, # added
       NeighborMesh,
       ReferenceDistance,
       closest_point_in_world,
       closest_point_in_body,
       separation_distance,
       simplex_penetration_distance,
       linear_combination,
       normal,  # PCA method
       normal_test  # PCA method unit test

include("tagged_points.jl")
include("simplices.jl")
# include("johnson_distance.jl")
include("neighbor_mesh.jl")
include("traits.jl")
include("gjk.jl")
include("reference_distance.jl")
include("johnson_distance.jl")
include("contact_normal.jl")

end # module
