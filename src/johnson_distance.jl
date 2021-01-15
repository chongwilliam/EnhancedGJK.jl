"""
    weights, wids, nvrtx = signed_volume(simplex, wids, nvrtx)

This function implements the Signed Volume distance sub-algorithm from
Montanari, Petrinijc: Improving the GJK Algorithm for faster and more
reliable distance queries between convex objects.  The advantage is to avoid
Johnson's sub-algorithm numerical instability in the Ax=b system for calculating
the next v_k vector.

M is number of points in simplex, and SVector{N,T} is size 3 and type Float64.
In the main algorithm, M is always fully populated (size = 4), and the weights
are set to zero for non-active simplex points. wids tracks the top nvrtx active vertices.
"""
function signed_volume(simplex::SVector{M,SVector{N,T}}, wids::SVector{M,P}, nvrtx::P) where {M,N,T,P}
    # Calculate weights and mininum active support set
    if nvrtx == 4
        weights, wids, nvrtx = S3D(simplex, wids)
        return weights, wids, nvrtx
    elseif nvrtx == 3
        weights, wids, nvrtx = S2D(simplex, wids)
        return weights, wids, nvrtx
    elseif nvrtx == 2
        weights, wids, nvrtx = S1D(simplex, wids)
        return weights, wids, nvrtx
    elseif nvrtx == 1
        weights = SVector{M,T}(1, 0, 0, 0)
        weights = order_weights(weights, wids)
        return weights, wids, nvrtx
    end
end

"""
S3D function for determining minimum active support set for tetrahedron and associated weights.
"""
function S3D(simplex::SVector{M,SVector{N,T}}, wids::SVector{M,P}) where {M,N,T,P}
    # Unpack simplex in wids order
    _s = view(simplex,1:4)

    # Calculate signed determinants to compare tetrahedron orientations
    M_41 = SMatrix{3,3,T}(_s[2][1],_s[2][2],_s[2][3],_s[3][1],_s[3][2],_s[3][3],_s[4][1],_s[4][2],_s[4][3])
    M_42 = SMatrix{3,3,T}(_s[1][1],_s[1][2],_s[1][3],_s[3][1],_s[3][2],_s[3][3],_s[4][1],_s[4][2],_s[4][3])
    M_43 = SMatrix{3,3,T}(_s[1][1],_s[1][2],_s[1][3],_s[2][1],_s[2][2],_s[2][3],_s[4][1],_s[4][2],_s[4][3])
    M_44 = SMatrix{3,3,T}(_s[1][1],_s[1][2],_s[1][3],_s[2][1],_s[2][2],_s[2][3],_s[3][1],_s[3][2],_s[3][3])
    C = SVector{4,T}(-det(M_41), det(M_42), -det(M_43), det(M_44))
    detM = sum(C)
    facets_test = compare_signs.(detM, C)

    # Calculate minimum support set and associated weights
    if facets_test[1] == facets_test[2] == facets_test[3] == facets_test[4] == true
        weights = C / detM
        new_wids = SVector{4,Int}(1, 2, 3, 4)
        return weights, new_wids, 4
    else
        d_min = Inf
        best_weights, best_wids, best_nvrtx = undef, undef, undef
        tmp_weights, tmp_wids, tmp_nvrtx = undef, undef, undef
        @inbounds for i = 1:4
            if facets_test[i] == false  # wrt s ordering
                if i == 1
                    new_wids = SVector{4,Int}(2, 3, 4, 1)
                elseif i == 2
                    new_wids = SVector{4,Int}(1, 3, 4, 2)
                elseif i == 3
                    new_wids = SVector{4,Int}(1, 2, 4, 3)
                elseif i == 4
                    new_wids = SVector{4,Int}(1, 2, 3, 4)
                end
                tmp_weights, tmp_wids, tmp_nvrtx = S2D(simplex, new_wids)  # weights are in order
                v = linear_combination(tmp_weights, simplex)
                v_len = dot(v, v)
                if v_len < d_min
                    best_weights, d_min, best_wids, best_nvrtx = tmp_weights, v_len, tmp_wids, tmp_nvrtx
                end
            end
        end
        return best_weights, best_wids, best_nvrtx
    end
end

"""
S2D function for determining minimum active support set for plane (triangle) and associated weights.
Modification to the original algorithm by rotating the plane to a flat orientation.
"""
function S2D(simplex::SVector{M,SVector{N,T}}, wids::SVector{M,P}) where {M,N,T,P}
    # Unpack indices
    _wids = view(wids,1:4)

    # Unpack simplex in wids order
    _s = view(simplex[wids],1:3)

    # Setup normals and origin projection
    AB = _s[2] - _s[1]
    AC = _s[3] - _s[1]
    n = cross(AB, AC)
    p = (dot(_s[1], n)*n) / dot(n, n)
    _p = view(p,1:3)

    # Calculated signed determinants
    # ABC
    xyABC = SMatrix{3,3,T}(_s[1][1], _s[1][2], 1, _s[2][1], _s[2][2], 1, _s[3][1], _s[3][2], 1)
    yzABC = SMatrix{3,3,T}(_s[1][2], _s[1][3], 1, _s[2][2], _s[2][3], 1, _s[3][2], _s[3][3], 1)
    xzABC = SMatrix{3,3,T}(_s[1][1], _s[1][3], 1, _s[2][1], _s[2][3], 1, _s[3][1], _s[3][3], 1)
    signed_area = SVector{3,T}(det(yzABC), det(xzABC), det(xyABC))
    _, I3 = findmax(abs.(signed_area))  # index to exclude
    if I3 == 1
        I1, I2 = 2, 3
    elseif I3 == 2
        I1, I2 = 1, 3
    elseif I3 == 3
        I1, I2 = 1, 2
    end
    # ABC
    nu_max = signed_area[I3]
    # ABP vs. ABC
    ABP = SMatrix{3,3,T}(_s[1][I1], _s[1][I2], 1, _s[2][I1], _s[2][I2], 1, _p[I1], _p[I2], 1)
    # BCP vs. BCA
    BCP = SMatrix{3,3,T}(_s[2][I1], _s[2][I2], 1, _s[3][I1], _s[3][I2], 1, _p[I1], _p[I2], 1)
    # ACP vs. ACB
    ACP = SMatrix{3,3,T}(_s[1][I1], _s[1][I2], 1, _s[3][I1], _s[3][I2], 1, _p[I1], _p[I2], 1)
    C = SVector{4,T}(det(BCP), -det(ACP), det(ABP), 0)  # wrt s ordering
    facets_test = compare_signs.(nu_max, C)

    # Calculate minimum support set and associated weights
    if facets_test[1] == facets_test[2] == facets_test[3]
        weights = C / nu_max
        weights = order_weights(weights, wids)
        return weights, wids, 3
    else
        d_min = Inf
        best_weights, best_wids, best_nvrtx = undef, undef, undef
        tmp_weights, tmp_wids, tmp_nvrtx = undef, undef, undef
        @inbounds for i = 1:3
            if facets_test[i] == false
                if i == 1
                    new_wids = SVector{4,P}(_wids[2], _wids[3], _wids[4], _wids[1])
                elseif i == 2
                    new_wids = SVector{4,P}(_wids[1], _wids[3], _wids[4], _wids[2])
                elseif i == 3
                    new_wids = SVector{4,P}(_wids[1], _wids[2], _wids[4], _wids[3])
                end
                tmp_weights, tmp_wids, tmp_nvrtx = S1D(simplex, new_wids)  # weights are in order
                v = linear_combination(tmp_weights, simplex)
                v_len = dot(v, v)
                if v_len < d_min
                    best_weights, d_min, best_wids, best_nvrtx = tmp_weights, v_len, tmp_wids, tmp_nvrtx
                end
            end
        end
        return best_weights, best_wids, best_nvrtx
    end
end

"""
S1D function for determining minimum active support set for line and associated weights.
Modification to the original algorithm with the calculation of the weights.
"""
function S1D(simplex::SVector{M,SVector{N,T}}, wids::SVector{M,P}) where {M,N,T,P}

    # Unpack indices
    _wids = view(wids,1:2)

    # Unpack simplex in wids order
    _s = view(simplex[_wids],:)

    # Calculate projection
    AB = _s[2] - _s[1]
    t = - dot(_s[1], AB) / dot(AB, AB)  # s(t) = (1-t)*s1 + t*s2

    if t < 0
        weights = SVector{4,T}(1, 0, 0, 0)
        weights = order_weights(weights, wids)
        return weights, wids, 1
    elseif t < 1
        weights = SVector{4,T}(1-t, t, 0, 0)
        weights = order_weights(weights, wids)
        return weights, wids, 2
    else
        weights = SVector{4,T}(0, 1, 0, 0)
        weights = order_weights(weights, wids)
        return weights, wids, 1
    end
end

"""
Helper function for sign comparison used in orientation determination.
"""
function compare_signs(a::T, b::T) where {T}
    if (a>0) && (b>0) || (a<0) && (b<0)
        return true
    else
        return false
    end
end

"""
Helper function to order weights based on wids.
"""
function order_weights(weights::SVector{M,T}, wids::SVector{M,P}) where {M,T,P}
    ordered_weights = MVector{M,T}(undef)
    ordered_weights[wids] = weights
    return SVector(ordered_weights)
end

num_johnson_subsets(simplex_length::Integer) = 2^simplex_length - 1

"""
Compute all subsets of the points in the simplex in a reliable order. The
order is arbitrary, but was chosen to match the implemenation in
S. Cameron, “Enhancing GJK: computing mininum and penetration distances between
convex polyhedra,”. Specifically, subset i contains point j iff the binary
representation of i has a one at the jth bit.
"""
function johnson_subsets(simplex_length::Integer)
    num_subsets = num_johnson_subsets(simplex_length)
    subsets = falses(simplex_length, num_subsets)
    for i in 1:num_subsets
        for j in 1:simplex_length
            subsets[j, i] = (i >> (j - 1)) & 1
        end
    end
    subsets
end

"""
    weights = projection_weights(simplex)

This function implements Johnson's distance subalgorithm, as described in
E. G. Gilbert, D. W. Johnson, and S. S. Keerthi, “A fast procedure for
computing the distance between complex objects in three-dimensional space,”
1988. Given a simplex (a length N+1 vector of points of dimension N), it
returns weights such that LinearAlgebra.dot(weights, simplex) yields the point in the convex
hull of the simplex which is closest to the origin.

This is the critical loop of the GJK algorithm, so it has been heavily optimized
to precompute, inline, and unroll as nuch of the algorithm as possible. For a
more readable (and nuch slower) implementation, see
projection_weights_reference()
"""
@generated function projection_weights(simplex::SVector{M, SVector{N, T}}; atol=sqrt(eps(T))) where {M,N,T}
    simplex_length = M
    num_subsets = num_johnson_subsets(M)
    subsets = johnson_subsets(M)
    complements = .!subsets

    expr = quote
        deltas = SVector(tuple($([zero(SVector{simplex_length, T}) for i = 1 : num_subsets]...)))
    end

    # Set the weight of every singleton subset to 1.
    for i in 1:simplex_length
        elements = [:(zero($T)) for j in 1:simplex_length]
        elements[i] = :(one($T))
        arg = :(SVector{$simplex_length, $T}(tuple($(elements...))))

        push!(expr.args, quote
            deltas = setindex(deltas, $(arg), $(2^(i - 1)))
        end)
    end

    for s in 1:(num_subsets - 1)
        k = first((1:M)[subsets[:,s]])
        push!(expr.args, quote
            viable = true
        end)
        for j in (1:M)[complements[:,s]]
            s2 = s + (1 << (j - 1))
            push!(expr.args, quote
                d = deltas[$s2][$j]
            end)
            for i in (1:M)[subsets[:,s]]
                push!(expr.args, quote
                    d += deltas[$s][$i] *
                    (LinearAlgebra.dot(simplex[$i], simplex[$k]) - LinearAlgebra.dot(simplex[$i], simplex[$j]))
                end)
            end
            push!(expr.args, quote
                deltas = setindex(deltas, setindex(deltas[$s2], d, $j), $s2)
            end)
            push!(expr.args, quote
                if d > atol # See #17
                    viable = false
                end
            end)
        end
        push!(expr.args, quote
            viable = viable && all($(Expr(:tuple,
            [:(deltas[$s][$i] >= 0) for i in (1:M)[subsets[:,s]]]...)))
            if viable
                return deltas[$s] ./ sum(deltas[$s])
            end
        end)
    end
    push!(expr.args, quote
        return deltas[$num_subsets] ./ sum(deltas[$num_subsets])
    end)
    return expr
end

function projection_weights_reference(simplex::SVector{M, SVector{N, T}}) where {M,N,T}
    num_subsets = 2^M - 1
    complements = falses(M, num_subsets)
    subsets = falses(M, num_subsets)
    for i in 1:num_subsets
        for j in 1:M
            subsets[j, i] = (i >> (j - 1)) & 1
            complements[j, i] = !subsets[j, i]
        end
    end

    deltas = zeros(MMatrix{M, num_subsets, T})
    viable = ones(MVector{num_subsets, Bool})

    for i in 1:M
        deltas[i, 2^(i - 1)] = 1
    end

    for s in 1:(num_subsets - 1)
        k = first((1:M)[subsets[:,s]])
        for j in (1:M)[complements[:,s]]
            s2 = s + (1 << (j - 1))
            delta_j_s2 = 0
            for i in (1:M)[subsets[:,s]]
                delta_j_s2 += deltas[i, s] * (dot(simplex[i], simplex[k]) - dot(simplex[i], simplex[j]))
            end
            if delta_j_s2 > 0
                viable[s] = false
            elseif delta_j_s2 < 0
                viable[s2] = false
            end
            deltas[j, s2] = delta_j_s2
        end
        if viable[s]
            return deltas[:,s] ./ sum(deltas[:,s])
        end
    end
    return deltas[:,end] ./ sum(deltas[:,end])
end
