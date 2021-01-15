struct Difference{PA, PB}
    a::PA
    b::PB
end

Base.:*(n::Number, d::Difference) = Difference(*(n, d.a), *(n, d.b))
Base.:+(d1::Difference, d2::Difference) = Difference(d1.a + d2.a, d1.b + d2.b)
any_inside(d::Difference) = Difference(any_inside(d.a), any_inside(d.b))

struct CollisionCache{GeomA, GeomB, M, D <: Difference}
    bodyA::GeomA
    bodyB::GeomB
    simplex_points::MVector{M, D}
end

function CollisionCache(geomA, geomB)
    N = dimension(typeof(geomA))
    @assert dimension(typeof(geomA)) == dimension(typeof(geomB))
    CollisionCache(N, geomA, geomB)
end

function CollisionCache(::Val{N}, geomA, geomB) where N
    interior_point = any_inside(Difference(geomA, geomB))
    simplex = MVector(ntuple(i -> interior_point, Val(N + 1)))
    CollisionCache(geomA, geomB, simplex)
end

function reset!(cache::CollisionCache)
    interior_point = any_inside(Difference(cache.bodyA, cache.bodyB))
    simplex = cache.simplex_points
    @inbounds for i in eachindex(simplex)
        simplex[i] = interior_point
    end
    cache
end

dimension(::Type{CollisionCache{G1, G2, M, D}}) where {G1,G2,M,D} = dimension(G1)

function support_vector_max(geometry::gt.GeometryPrimitive, direction, initial_guess::Tagged)
    best_pt, score = gt.support_vector_max(geometry, gt.Vec(direction))
    Tagged(SVector(best_pt))
end

function support_vector_max(pt::gt.Vec{N, T}, direction, initial_guess::Tagged) where {N,T}
    Tagged(SVector(pt))
end

function support_vector_max(pt::gt.Point{N, T}, direction, initial_guess::Tagged) where {N,T}
    Tagged(SVector(pt))
end

function support_vector_max(simplex::Union{gt.AbstractSimplex, gt.AbstractFlexibleGeometry}, direction, initial_guess::Tagged)
    best_pt, score = gt.support_vector_max(simplex, gt.Vec(direction))
    Tagged(SVector(best_pt))
end

function support_vector_max(mesh::gt.HomogenousMesh{gt.Point{N, T}}, direction, initial_guess::Tagged) where {N,T}
    best_arg, best_value = gt.argmax(x-> dot(SVector(x), direction), gt.vertices(mesh))
    best_vec = SVector(best_arg)
    Tagged(best_vec)
end

any_inside(pt::SVector) = Tagged(pt)
support_vector_max(pt::SVector, direction, initial_guess::Tagged) = Tagged(pt)

function transform_simplex(cache::CollisionCache, poseA, poseB)
    transform_simplex(dimension(typeof(cache)), cache, poseA, poseB)
end

@generated function transform_simplex(::Val{N}, cache::CollisionCache, poseA, poseB) where {N}
    transform_simplex_impl(N, cache, poseA, poseB)
end

function transform_simplex_impl(N, cache, poseA, poseB)
    Expr(:call, :(SVector),
        [:((poseA(value(cache.simplex_points[$i].a)) -
            poseB(value(cache.simplex_points[$i].b)))) for i in 1:(N + 1)]...)
end

# Note: it looks like this can be replaced with transpose(weights) * points in Julia 1.3 (before that, it's a lot slower)
# @generated function linear_combination(weights::StaticVector{N}, points::StaticVector{N}) where {N}
#     expr = :(weights[1] * points[1])
#     for i = 2 : N
#         expr = :($expr + weights[$i] * points[$i])
#     end
#     return quote
#         Base.@_inline_meta
#         $expr
#     end
# end

# Updated for Julia 1.5
function linear_combination(weights::StaticVector{N}, points::StaticVector{N}) where {N}
    return weights' * points
end

struct GJKResult{M, N, T}
    simplex::SVector{M, SVector{N, T}}
    weights::SVector{M, T}
    in_collision::Bool
    closest_point_in_body::Difference{SVector{N, T}, SVector{N, T}}
    nvrtx::Int
    wids::SVector{M, Int}
    iterations::Int
    termination::Int
end

function Base.getproperty(result::GJKResult, sym::Symbol)
    if sym === :signed_distance
        @warn "The `signed_distance` field was removed. Please use the `separation_distance` and `simplex_penetration_distance` functions."
        if result.in_collision
            return -simplex_penetration_distance(result)
        else
            return separation_distance(result)
        end
    else
        return Core.getfield(result, sym)
    end
end

""" Calculate active simplex points of privileged object for neighbor matching in implicit surface in body frame.  """
function body_simplex(cache::CollisionCache, wids::SVector{M, Int}, nvrtx::Int) where {M}
    simplex_points = [value(cache.simplex_points[i].a) for i in wids[1:nvrtx]]
    # T = typeof(simplex_points[1][2])  # Number type
    return SVector{nvrtx, SVector{3, Float64}}(simplex_points)
end

closest_point_in_world(result::GJKResult) = linear_combination(result.simplex, result.weights)
closest_point_in_body(result::GJKResult) = result.closest_point_in_body
separation_distance(result::GJKResult) = (@assert !result.in_collision; norm(closest_point_in_world(result)))
simplex_penetration_distance(result::GJKResult) = penetration_distance(result.simplex)

"""
Algorithm referenced from Bergen "Collision Detection in Interactive 3D Environments" textbook, page 145.
"""
function gjk!(cache::CollisionCache,
              poseA::Transformation,
              poseB::Transformation,
              max_iter=10,
              eps_tol=1e-16,  # improvement point tolerance
              eps_rel=1e-16)  # collision tolerance

    rotAinv = transform_deriv(inv(poseA), 0)
    rotBinv = transform_deriv(inv(poseB), 0)
    simplex = transform_simplex(cache, poseA, poseB)
    iter = 1

    # Initialize containers
    weights = SVector{4,Float64}(1,0,0,0)
    wids = SVector{4,Int}(1,2,3,4)  # track active simplex points returned from signed_volume (top nvrtx are active)
    prev_wids = wids  # tracks repetition termination condition
    nvrtx = 0  # number of active vertices in simplex
    best_point = simplex[1]  # initialize to first vertex

    while true
        # Compute the support function to find w_k and replace in simplex
        direction = -best_point
        direction_in_A = rotAinv * direction
        direction_in_B = rotBinv * direction

        starting_vertex_index = 1
        starting_vertex = cache.simplex_points[starting_vertex_index]
        starting_vertex_score =
            dot(value(starting_vertex.a), direction_in_A) +
            dot(value(starting_vertex.b), direction_in_B)
        for j in 2:length(cache.simplex_points)
            candidate = cache.simplex_points[j]
            candidate_score =
                dot(value(candidate.a), direction_in_A) +
                dot(value(candidate.b), direction_in_B)
            if candidate_score > starting_vertex_score
                starting_vertex_score = candidate_score
                starting_vertex = candidate
            end
        end

        improved_vertex = Difference(
            support_vector_max(cache.bodyA, direction_in_A, starting_vertex.a),
            support_vector_max(cache.bodyB, -direction_in_B, starting_vertex.b))
        improved_point = poseA(value(improved_vertex.a)) - poseB(value(improved_vertex.b))
        # score = dot(improved_point, direction)  # best_point = v_k; improved_point = w_k

        # println("***")
        # println("pre-update improved point: ", improved_point)

        # Primary termination conditions
        v_len = dot(best_point, best_point)
        if iter == max_iter
            # Max iterations
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, false, closest_point_in_body, nvrtx, wids, iter, 1)
        elseif v_len - dot(best_point, improved_point) <= eps_rel*v_len
            # w_k is sufficiently close enough to v_k to not improve v
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, false, closest_point_in_body, nvrtx, wids, iter, 2)
        elseif improved_point in simplex[prev_wids]
            # w_k ∈ W_(k-1) U w_(k-1)
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, false, closest_point_in_body, nvrtx, wids, iter, 3)
        else
            # Add improved point to the simplex set
            nvrtx += 1
            index_to_replace = wids[nvrtx]
            cache.simplex_points[index_to_replace] = improved_vertex
            simplex = setindex(simplex, improved_point, index_to_replace)
            prev_wids = wids[1:nvrtx]
            # println("pre-update wids: ", wids)
            # println("pre-update vrtx: ", nvrtx)
            # println("pre-update previous wids: ", prev_wids)
        end

        # Calculate new weights and support set W
        weights, wids, nvrtx = signed_volume(simplex, wids, nvrtx)
        # min_weight, index_to_replace = findmin(weights)

        # Secondary termination conditions
        # best_point = linear_combination(weights, simplex)
        # w_set = simplex[wids[1:nvrtx]]
        # w_max = maximum(dot(w_set, w_set))
        #
        # if nvrtx == 4 || dot(best_point, best_point) <= eps_tol*w_max
        #     closest_point_in_body = linear_combination(weights, cache.simplex_points)
        #     return GJKResult(simplex, weights, true, closest_point_in_body, nvrtx, wids, iter, 4)
        # end

        best_point = linear_combination(weights, simplex)
        if nvrtx == 4
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, true, closest_point_in_body, nvrtx, wids, iter, 4)
        end

        # println("post-update nvrtx: ", nvrtx)
        # println("post-update weights: ", weights)
        # println("post-update wids: ", wids)

        iter += 1
    end
end

function penetration_distance(simplex)
    best_dist_squared = nothing
    for i in eachindex(simplex)
        face = simplex_face(simplex, i)
        weights = projection_weights(face)
        closest_point = linear_combination(weights, face)
        dist_squared = closest_point ⋅ closest_point
        if best_dist_squared === nothing || dist_squared < best_dist_squared
            best_dist_squared = dist_squared
        end
    end
    return sqrt(best_dist_squared)
end

function gjk(geomA, geomB,
             poseA::Transformation=IdentityTransformation(),
             poseB::Transformation=IdentityTransformation())
    cache = CollisionCache(geomA, geomB)
    gjk!(cache, poseA, poseB)
end

""" Using ForwardDiff for gradients """
function gjk_diff(q::AbstractVector{T}) where {T}  # q is generalized coordinates
    RotA = composed(RotZ(q[6]), composed(RotY(q[5]), RotX(q[4])))
    TransA = Translation(q[1:3])
    # RotB = RotZ(RotY(RotX(q[4:6])))
    # TransB = Translation(q[1:3])
    poseA = composed(RotA, TransA)
    poseB = IdentityTransformation
    result = gjk(c1, c2, poseA, poseB)
    return separation_distance(result)
end

#---
### Original functions ###
function gjk_original!(cache::CollisionCache,
              poseA::Transformation,
              poseB::Transformation,
              max_iter=100,
              atol=1e-6)
    rotAinv = transform_deriv(inv(poseA), 0)
    rotBinv = transform_deriv(inv(poseB), 0)
    simplex = transform_simplex(cache, poseA, poseB)
    iter = 1

    while true
        weights = projection_weights(simplex)
        min_weight, index_to_replace = findmin(weights)
        best_point = linear_combination(weights, simplex)
        if min_weight > 0
            # Nominally in collision; but check for numerical issues
            separation_squared = best_point ⋅ best_point
            @assert separation_squared < sqrt(1_000_000 * eps(typeof(separation_squared)))
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            # return GJKResult(simplex, weights, true, closest_point_in_body, iter, simplex, poseA, poseB)
            return GJKResult(simplex, weights, true, closest_point_in_body, 0, SVector{4,Int}(zeros(4,1)), iter, 0)
        end

        direction = -best_point
        direction_in_A = rotAinv * direction
        direction_in_B = rotBinv * direction

        starting_vertex_index = 1
        starting_vertex = cache.simplex_points[starting_vertex_index]
        starting_vertex_score =
            dot(value(starting_vertex.a), direction_in_A) +
            dot(value(starting_vertex.b), direction_in_B)
        for j in 2:length(cache.simplex_points)
            candidate = cache.simplex_points[j]
            candidate_score =
                dot(value(candidate.a), direction_in_A) +
                dot(value(candidate.b), direction_in_B)
            if candidate_score > starting_vertex_score
                starting_vertex_score = candidate_score
                starting_vertex = candidate
            end
        end

        improved_vertex = Difference(
            support_vector_max(cache.bodyA, direction_in_A, starting_vertex.a),
            support_vector_max(cache.bodyB, -direction_in_B, starting_vertex.b))
        improved_point = poseA(value(improved_vertex.a)) - poseB(value(improved_vertex.b))
        score = dot(improved_point, direction)
        if score <= dot(best_point, direction) + atol || iter >= max_iter
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, false, closest_point_in_body, 0, SVector{4,Int}(zeros(4,1)), iter, 0)
        else
            cache.simplex_points[index_to_replace] = improved_vertex
            simplex = setindex(simplex, improved_point, index_to_replace)
            # simplex[index_to_replace] = improved_point
        end

        iter += 1
    end
end

function gjk_original(geomA, geomB,
             poseA::Transformation=IdentityTransformation(),
             poseB::Transformation=IdentityTransformation())
    cache = CollisionCache(geomA, geomB)
    gjk_original!(cache, poseA, poseB)
end
