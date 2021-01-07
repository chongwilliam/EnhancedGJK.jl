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
    weights, wids = signed_volume(simplex, weights, wids, nvrtx)

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
        S3D(simplex, weights, wids, nvrtx)
        return weights, wids, nvrtx
    elseif nvrtx == 3
        S2D(simplex, weights, wids, nvrtx)
        return weights, wids, nvrtx
    elseif nvrtx == 2
        S1D(simplex, weights, wids, nvrtx)
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
    B = MVector{4,T}(B1, B2, B3, B4)
    detM = sum(B)
    FacetsTest = CompareSigns.(detM, B)

    # Calculate minimum support set and associated weights
    if FacetsTest .== true
        weights = B ./ detM
        nvrtx = 4
        return weights, wids, nvrtx
    else
        d_min = Inf
        best_weights, best_wids = MVector{4,T}(zeros(4,1)), MVector{4,P}(zeros(4,1))
        for i = 1:4
            if FacetsTest[i] == false
                nvrtx = 3
                if i == 1
                    wids[1], wids[2], wids[3], wids[4] = s2, s3, s4, s1
                elseif i == 2
                    wids[1], wids[2], wids[3], wids[4] = s1, s3, s4, s2
                elseif i == 3
                    wids[1], wids[2], wids[3], wids[4] = s1, s2, s4, s3
                end
                S2D!(simplex, weights, wids, nvrtx)  # exclude index i
                v = linear_combination(weights, simplex)
                d = dot(v, v)
                if d < d_min
                    best_weights = weights
                    d_min = d
                    best_wids = wids
                end
            end
        end
        weights, wids, nvrtx = best_weights, best_wids, nvrtx  # get best values
    end
end

"""
S2D function for determining minimum active support set for plane (triangle) and associated weights.
Modification to the original algorithm by rotating the plane to a flat orientation.
"""
function S2D!(simplex::SVector{M,SVector{N,T}}, weights::MVector{M,T}, wids::MVector{M,P}, nvrtx::P) where {M,N,T,P}
    # Unpack indices
    s1, s2, s3 = wids[1:3]

    # Setup normals and origin projection
    a, b, c = s1, s2, s3
    s21 = b - a
    s31 = c - a
    n = LinearAlgebra.cross(s21, s31)
    p0 = (dot(a, n))*n / dot(n, n)

    # Method using DCM rotation matrix
    t1_hat = s21 / norm(s21)
    n_hat = n / norm(n)
    t2_hat = cross(n_hat, t1_hat)

    R = SMatrix{3,3,T}(t1_hat[1], t1_hat[2], t1_hat[3], t2_hat[1], t2_hat[2], t2_hat[3], n_hat[1], n_hat[2], n_hat[3])  # rotate from plane to N
    a = DCM' * a
    b = DCM' * b
    c = DCM' * c
    p0 = DCM' * p0
    nu_max = a[1]*b[2] - a[2]*b[1] - a[1]*c[2] + a[2]*c[1] + b[1]*c[2] - b[2]*c[1]  # signed area of triangle

    # Discard the j-th coordinate to compare orientations along the plane tangent directions
    # a_proj = SVector{2,T}(a[1:end .!= 3])
    # b_proj = SVector{2,T}(b[1:end .!= 3])
    # c_proj = SVector{2,T}(c[1:end .!= 3])
    # p_proj = p0[1:end .!= 3]

    # Compute signed determinants
    # B = MVector{3,T}(0, 0, 0)
    # s_mat = SMatrix{3,3,T}(p_proj[1], p_proj[2], 1, b_proj[1], b_proj[2], 1, c_proj[1], c_proj[2], 1)  # determinant of this matrix
    # B1 = b_proj[1]*c_proj[2] - b_proj[2]*c_proj[1] - b_proj[1]*p_proj[2] + b_proj[2]*p_proj[1] + c_proj[1]*p_proj[2] - c_proj[2]*p_proj[1]
    B1 = b[1]*c[2] - b[2]*c[1] - b[1]*p[2] + b[2]*p[1] + c[1]*p[2] - c[2]*p[1]

    # s_mat = SMatrix{3,3,T}(a_proj[1], a_proj[2], 1, p_proj[1], p_proj[2], 1, c_proj[1], c_proj[2], 1)  # determinant of this matrix
    # B2 = a_proj[2]*c_proj[1] - a_proj[1]*c_proj[2] + a_proj[1]*p_proj[2] - a_proj[2]*p_proj[1] - c_proj[1]*p_proj[2] + c_proj[2]*p_proj[1]
    B2 = a[2]*c[1] - a[1]*c[2] + a[1]*p[2] - a[2]*p[1] - c[1]*p[2] + c[2]*p[1]

    # s_mat = SMatrix{3,3,T}(a_proj[1], a_proj[2], 1, b_proj[1], b_proj[2], 1, p_proj[1], p_proj[2], 1)  # determinant of this matrix
    # B3 = a_proj[1]*b_proj[2] - a_proj[2]*b_proj[1] - a_proj[1]*p_proj[2] + a_proj[2]*p_proj[1] + b_proj[1]*p_proj[2] - b_proj[2]*p_proj[1]
    B3 = a[1]*b[2] - a[2]*b[1] - a[1]*p[2] + a[2]*p[1] + b[1]*p[2] - b[2]*p[1]

    B = MVector{3,T}(B1, B2, B3)
    FacetsTest = CompareSigns.(nu_max, B)

    # FacetsTest = [0, 0, 0]
    # for i = 1:3
        # FacetsTest[i] = CompareSigns(nu_max, B[i])
    # end

    # Calculate minimum support set and associated weights
    if FacetsTest .== true
        nvrtx = 3
        weights[s1], weights[s2], weights[s3] = B ./ nu_max
        weights[s4] = 0
    else
        # print("S2D Condition 2", "\n")
        d = 1e5
        for i = 1:3
            if FacetsTest[i] == 0
                new_ind = MVector(ind)  # bit array carried over from function input
                new_ind[actual_ind[i]] = 0  # remove the actual i-th coordinate
                lam_new, w_new = S1D(simplex, new_ind)
                v_new = linear_combination(lam_new, simplex)
                d_new = LinearAlgebra.dot(v_new, v_new)
                if d_new < d
                    lam = lam_new
                    d = d_new
                    w_ind = w_new
                end
            end
        end
        return lam, w_ind
    end
end
"""
S1D function for determining minimum active support set for line and associated weights.
Modification to the original algorithm with the calculation of the weights.
"""
function S1D(simplex::SVector{M, SVector{N, T}}, ind::MVector{M,P}) where {M,N,T,P}
    # Unpack with indices
    actual_ind = MVector{2,P}(0, 0)  # stores actual indices for the vertices/weights
    s = []  # stores simplices used in the projection
    cnt = 1
    for i = 1:4
        if ind[i] == 1
            push!(s, simplex[i])
            actual_ind[cnt] = i
            cnt += 1
        end
    end

    # Project origin onto line segment AB
    a, b = s[1], s[2]
    t = b - a
    pt = b + (LinearAlgebra.dot(b, t)/LinearAlgebra.dot(t, t)) * (-t)

    # Determine safest projection axis
    nu_max, I_ind = 0, 1
    for i = 1:3
        nu = t[i]
        if abs(nu) > abs(nu_max)
            nu_max = nu
            I_ind = i
        end
    end

    a_proj, b_proj, pt_proj = a[I_ind], b[I_ind], pt[I_ind]
    B = MVector{2,T}(0, 0)
    B[1] = pt_proj - a_proj
    B[2] = b_proj - pt_proj

    FacetsTest = [0, 0]
    FacetsTest[1] = CompareSigns(nu_max, B[1])
    FacetsTest[2] = CompareSigns(nu_max, B[2])

    # Calculate minimum support set and associated weights
    lam = MVector{4,T}(0, 0, 0, 0)
    w_ind = MVector{4,P}(0, 0, 0, 0)
    if (FacetsTest[1] == 1 && FacetsTest[2] == 1)
        # Origin is between A and B
        # print("S1D Condition 1", "\n")
        lam[actual_ind[1]] = B[2] / nu_max
        lam[actual_ind[2]] = B[1] / nu_max
        w_ind[actual_ind[1]] = 1
        w_ind[actual_ind[2]] = 1
        return lam, w_ind
    else
        # Select vertex A or B based on distance from origin
        if LinearAlgebra.dot(a, a) < LinearAlgebra.dot(b, b)
            # print("S1D Condition 2", "\n")
            lam[actual_ind[1]] = 1
            w_ind[actual_ind[1]] = 1
            return lam, w_ind
        else
            # print("S1D Condition 3", "\n")
            lam[actual_ind[2]] = 1
            w_ind[actual_ind[2]] = 1
            return lam, w_ind
        end
    end
end

"""
Helper function for sign comparison used in orientation determination.
"""
function CompareSigns(a::T, b::T) where {T}
    return (a>0) && (b>0)
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
