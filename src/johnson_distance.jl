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
    weights, wids, nvrtx = signed_volume(simplex, weights, wids, nvrtx)

This function implements the Signed Volume distance sub-algorithm from
Montanari, Petrinijc: Improving the GJK Algorithm for faster and more
reliable distance queries between convex objects.  The advantage is to avoid
Johnson's sub-algorithm numerical instability in the Ax=b system for calculating
the next v_k vector.

M is number of points in simplex, and SVector{N,T} is size 3 and type Float64.
In the main algorithm, M is always fully populated (size = 4), and the weights
are set to zero for non-active simplex points.
"""
function signed_volume(simplex::SVector{M,SVector{N,T}}, weights::MVector{M,T}, wids::MVector{M,P}, nvrtx::P) where {M,N,T,P}
    # Calculate weights and mininum active support set
    if nvrtx == 4
        weights, wids, nvrtx = S3D(simplex, weights, wids, nvrtx)
        # println("S3D: ", nvrtx, weights)
        return weights, wids, nvrtx
    elseif nvrtx == 3
        weights, wids, nvrtx = S2D(simplex, weights, wids, nvrtx)
        # println("S2D: ", nvrtx, weights)
        return weights, wids, nvrtx
    elseif nvrtx == 2
        weights, wids, nvrtx = S1D(simplex, weights, wids, nvrtx)
        # println("S1D: ", nvrtx, weights)
        return weights, wids, nvrtx
    elseif nvrtx == 1
        return weights, wids, nvrtx
    end
end

"""
S3D function for determining minimum active support set for tetrahedron and associated weights.
"""
function S3D(simplex::SVector{M,SVector{N,T}}, weights::MVector{M,T}, wids::MVector{M,P}, nvrtx::P) where {M,N,T,P}
    # Unpack indices
    s1, s2, s3, s4 = wids

    # Calculate signed determinants to compare tetrahedron orientations
    a, b, c, d = simplex[s1], simplex[s2], simplex[s3], simplex[s4]
    B1 = -1*(b[1]*c[2]*d[3] - b[1]*c[3]*d[2] - b[2]*c[1]*d[3] + b[2]*c[3]*d[1] + b[3]*c[1]*d[2] - b[3]*c[2]*d[1])
    B2 = a[1]*c[2]*d[3] - a[1]*c[3]*d[2] - a[2]*c[1]*d[3] + a[2]*c[3]*d[1] + a[3]*c[1]*d[2] - a[3]*c[2]*d[1]
    B3 = -1*(a[1]*b[2]*d[3] - a[1]*b[3]*d[2] - a[2]*b[1]*d[3] + a[2]*b[3]*d[1] + a[3]*b[1]*d[2] - a[3]*b[2]*d[1])
    B4 = a[1]*b[2]*c[3] - a[1]*b[3]*c[2] - a[2]*b[1]*c[3] + a[2]*b[3]*c[1] + a[3]*b[1]*c[2] - a[3]*b[2]*c[1]
    B = SVector{4,T}(B1, B2, B3, B4)
    detM = sum(B)
    FacetsTest = CompareSigns.(detM, B)

    # Calculate minimum support set and associated weights
    if FacetsTest[1] == FacetsTest[2] == FacetsTest[3] == FacetsTest[4] == true
        weights[s1], weights[s2], weights[s3], weights[s4] = B[1]/detM, B[2]/detM, B[3]/detM, B[4]/detM
        return weights, wids, 4
    else
        d_min = Inf
        best_weights, best_wids, best_nvrtx = MVector{4,T}(zeros(4,1)), MVector{4,P}(zeros(4,1)), Int(0)
        for i = 1:4
            if FacetsTest[i] == false
                nvrtx = 3
                if i == 1
                    wids[1], wids[2], wids[3], wids[4] = s2, s3, s4, s1
                elseif i == 2
                    wids[1], wids[2], wids[3], wids[4] = s1, s3, s4, s2
                elseif i == 3
                    wids[1], wids[2], wids[3], wids[4] = s1, s2, s4, s3
                elseif i == 4
                    wids[1], wids[2], wids[3], wids[4] = s1, s2, s3, s4
                end
                S2D(simplex, weights, wids, nvrtx)  # exclude index i
                v = linear_combination(weights, simplex)
                d = dot(v, v)
                if d < d_min
                    best_weights = weights
                    d_min = d
                    best_wids = wids
                    best_nvrtx = nvrtx
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
function S2D(simplex::SVector{M,SVector{N,T}}, weights::MVector{M,T}, wids::MVector{M,P}, nvrtx::P) where {M,N,T,P}
    # Unpack indices
    s1, s2, s3, s4 = wids

    # Setup normals and origin projection
    a, b, c = simplex[s1], simplex[s2], simplex[s3]
    s21 = b - a
    s31 = c - a
    n = cross(s21, s31)
    p = (dot(a, n))*n / dot(n, n)

    # Method using DCM rotation matrix
    t1_hat = s21 / norm(s21)
    n_hat = n / norm(n)
    t2_hat = cross(n_hat, t1_hat)

    R = SMatrix{3,3,T}(t1_hat[1], t1_hat[2], t1_hat[3], t2_hat[1], t2_hat[2], t2_hat[3], n_hat[1], n_hat[2], n_hat[3])  # rotate from plane to N
    a, b, c, p = R'*a, R'*b, R'*c, R'*p
    nu_max = a[1]*b[2] - a[2]*b[1] - a[1]*c[2] + a[2]*c[1] + b[1]*c[2] - b[2]*c[1]  # signed area of triangle

    # Compute signed determinants
    B1 = b[1]*c[2] - b[2]*c[1] - b[1]*p[2] + b[2]*p[1] + c[1]*p[2] - c[2]*p[1]
    B2 = a[2]*c[1] - a[1]*c[2] + a[1]*p[2] - a[2]*p[1] - c[1]*p[2] + c[2]*p[1]
    B3 = a[1]*b[2] - a[2]*b[1] - a[1]*p[2] + a[2]*p[1] + b[1]*p[2] - b[2]*p[1]
    B = SVector{3,T}(B1, B2, B3)
    FacetsTest = CompareSigns.(nu_max, B)

    # Calculate minimum support set and associated weights
    if FacetsTest[1] == FacetsTest[2] == FacetsTest[3] == true
        weights[s1], weights[s2], weights[s3] = B[1]/nu_max, B[2]/nu_max, B[3]/nu_max
        weights[s4] = 0
        return weights, wids, 3
    else
        d_min = Inf
        best_weights, best_wids, best_nvrtx = MVector{4,T}(zeros(4,1)), MVector{4,P}(zeros(4,1)), Int(0)
        for i = 1:3
            if FacetsTest[i] == false
                nvrtx = 2
                if i == 1
                    wids[1], wids[2], wids[3] = s2, s3, s1
                elseif i == 2
                    wids[1], wids[2], wids[3] = s1, s3, s2
                elseif i == 3
                    wids[1], wids[2], wids[3] = s1, s2, s3
                end
                S1D(simplex, weights, wids, nvrtx)
                v = linear_combination(weights, simplex)
                d = dot(v, v)
                if d < d_min
                    best_weights = weights
                    d_min = d
                    best_wids = wids
                    best_nvrtx = nvrtx
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
function S1D(simplex::SVector{M,SVector{N,T}}, weights::MVector{M,T}, wids::MVector{M,P}, nvrtx::P) where {M,N,T,P}
    # Unpack with indices
    s1, s2, s3, s4 = wids

    # Project origin onto line segment AB
    a, b = simplex[s1], simplex[s2]
    t = b - a
    pt = b + (dot(b, t)/dot(t, t)) * (-t)

    # Determine safest projection axis
    nu_max, ind = findmax(abs.(t))
    B = SVector{2,T}(pt[ind] - a[ind], b[ind] - pt[ind])
    FacetsTest = CompareSigns.(nu_max, B)

    # Calculate minimum support set and associated weights
    if FacetsTest[1] == FacetsTest[2] == true
        weights[s1], weights[s2] = B[2]/nu_max, B[1]/nu_max
        weights[s3], weights[s4] = 0, 0
        return weights, wids, 2
    else
        # Select vertex A or B based on distance from origin
        if dot(a, a) < dot(b, b)
            weights .= 0
            weights[s1] = 1
            return weights, wids, 1
        else
            weights .= 0
            weights[s2] = 1
            wids[1], wids[2] = s2, s1
            return weights, wids, 1
        end
    end
end

"""
Helper function for sign comparison used in orientation determination.
"""
function CompareSigns(a::T, b::T) where {T}
    if (a>0) && (b>0) || (a<0) && (b<0)
        return true
    else
        return false
    end
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
