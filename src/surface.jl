""" Construct variational implicit surface given unordered set of points in R3 """

# struct surface{M,N,T}
#     centers::MVector{M,SVector{N,T}}
#     distances::MVector{M,T}
#     coeff::SVector{M,T}
#     A::MMatrix{P,P,T}
#     b::MVector{P,T}
# end

""" Duchon's interpolant """
function rbf(x, y)
    return norm(x-y)^3
end

# """ Form interpolation matrix A """
# function interpolation_matrix!(A::MMatrix{P,P,T}, centers::SVector{M,SVector{N,T}}) where {M,N,P,T}
#     for i = 1:M
#         A[i,P-2:P] = centers[i]
#         A[P-2:P,i] = centers[i]
#         for j = 1:M
#             A[i,j] = rbf(centers[i],centers[j])
#         end
#     end
#     A[M+1,1:M] .= 1
#     A[1:M,M+1] .= 1
#     A[M+1:P,M+1:P] .= 0
# end

""" Form interpolation matrix A """
function interpolation_matrix!(A::MMatrix{P,P,T}, centers::SVector{M,SVector{N,T}}) where {M,N,P,T}
    @inbounds for i = 1:M
        # Upper right block
        A[i,P-2] = centers[i][1]
        A[i,P-1] = centers[i][2]
        A[i,P] = centers[i][3]
        # Lower left block
        A[P-2,i] = centers[i][1]
        A[P-1,i] = centers[i][2]
        A[P,i] = centers[i][3]
        # Upper left block
        @inbounds for j = 1:M
            A[i,j] = rbf(centers[i],centers[j])
        end
        # Column of ones
        A[M+1,i] = 1
        #  Row of ones
        A[i,M+1] = 1
    end
    # Bottom right 4x4 block of zeros
    @inbounds for i = 1:4
        @inbounds for j = 1:4
            A[M+i,M+j] = 0
        end
    end
    return nothing
end

""" Form b vector corresponding to signed distance values """
function fit_vector!(b::MVector{M,T}, distances::SVector{N,T}) where {M,N,T}
    @inbounds for i = 1:N
        b[i] = distances[i]
    end
    @inbounds for i = 1:4
        b[N+i] = 0
    end
    return nothing
end

""" Evaluate implicit surface value given point x (single evaluation) """
function eval_surface(coeff::SVector{M,T}, centers::SVector{N,SVector{P,T}}, x) where {M,N,P,T}
    # f = c[M-3] + c[M-2]*x[1] + c[M-1]*x[2] + c[M]*x[3]
    f = coeff[M-3]
    f += dot(SVector{3,T}(coeff[M-2], coeff[M-1], coeff[M]), x)
    @inbounds for i = 1:N
        f += coeff[i]*rbf(x, centers[i])
    end
    return f
end

""" Main call """
function implicit_surface(verts::Vector{gt.Point{3,T}},
                          points::SVector{N,SVector{3,T}},
                          signed_distances::SVector{N,T}) where {N,T}
    centers, distances = construct_centers(verts, points, signed_distances)
    coeff = construct_surface(centers, distances)
    return coeff, centers
end

""" Construct implicit surface """
function construct_surface(centers::SVector{M,SVector{N,T}}, distances::SVector{M,T}) where {M,N,T}
    A = MMatrix{M+4,M+4,T}(undef)
    b = MVector{M+4,T}(undef)
    interpolation_matrix!(A, centers)
    fit_vector!(b, distances)
    coeff = A \ b
    return SVector{M+4,Float64}(coeff)
end

""" Construct centers of implicit surface """
function construct_centers(verts::Vector{gt.Point{3,T}},
                           points::SVector{N,SVector{3,T}},
                           point_distances::SVector{N,T}) where {N,T}
    M = length(verts)
    centers = MVector{M+N,SVector{3,Float64}}(undef)
    distances = MVector{M+N,Float64}(undef)
    @inbounds for i = 1:M
        centers[i] = SVector{3,Float64}(verts[i])
        distances[i] = 0
    end
    @inbounds for i = 1:N
        centers[M+i] = points[i]
        distances[M+i] = point_distances[i]
    end
    return SVector(centers), SVector(distances)
end

# """ Create initial rbf centers from geometry, where N is the number of environment points """
# function init_centers(verts::Vector{GeometryTypes.Point{3,T}}, N::Int) where {T}
#     M = length(verts)
#     centers = MVector{M+N,SVector{3,Float64}}(undef)
#     distances = MVector{M+N,Float64}(undef)
#     for i = 1:M
#         centers[i] = SVector{3,Float64}(verts[i])
#         distances[i] = 0
#     end
#     return centers, distances
# end

# """ Add points and distances to centers """
# function add_centers!(centers::MVector{M,SVector{3,T}}, distances::MVector{M,T}, points::SVector{N,SVector{3,T}}, point_distances::SVector{N,T}) where {M,N,T}
#
#     add_centers = MVector{M+N, SVector{3,Float64}}(undef)
#     add_dist = MVector{M+N, Float64}(zeros(M+N))
#     for i = 1:N
#         centers[i] = SVector{3,Float64}(verts[i])
#     end
#     for i = 1:M
#         centers[N+i] = points[i]
#         dist[N+i] = distances[i]
#     end
#     return SVector(centers), SVector(dist)
# end

# Functions for surface fitting
""" Get ∇f(x,y,z) and H(x,y,z) at given eval point (x,y,z) """
function surface_gradient!(∇f,
                           coeff::SVector{M,T},
                           centers::SVector{N,SVector{P,T}},
                           point) where {M,N,P,T}
    ForwardDiff.gradient!(∇f, x -> eval_surface(coeff, centers, x), point)
    # return ∇f
    # H = ForwardDiff.hessian(x -> eval_surface(coeff, centers, x), point)
    # return ∇f, H
    # result = ForwardDiff.gradient!(result, x -> eval_surface(coeff, centers, x), point)
    # result = ForwardDiff.hessian!(result, x -> eval_surface(coeff, centers, x), point)
    return nothing
end

""" First order newton method to project grid points onto zero level set """
function fzero(coeff, centers, grid, ∇f; max_iter=20, atol=1e-6, xtol=1e-3)
    n_grid = length(grid)
    zero_points = MVector{n_grid, SVector{3,Float64}}(undef)
    cnt = 1
    for i = 1:n_grid
        iter = 1
        x = grid[i]
        while true
            y = eval_surface(coeff, centers, x)
            if iter == max_iter
                break
            elseif dot(y, y) < atol
                zero_points[cnt] = SVector{3,Float64}(x)
                cnt += 1
                break
            else
                surface_gradient!(∇f, coeff, centers, x)
                Δx = y * (∇f / dot(∇f, ∇f))
                println(Δx)
                if dot(Δx, Δx) < xtol
                    break
                else
                    x -= Δx
                end
                iter += 1
            end
        end
    end
    return zero_points, cnt-1
end

# Meshing functions
""" Return points and faces of the isosurface of the implicit surface """
function make_geom(coeff::SVector{M,T}, centers::SVector{N, SVector{P,T}}, origin::SVector{3,T}, dims::SVector{3,T}, pts::Int) where {M,N,P,T}
    function f_eval(x)
        f = coeff[N-3] + coeff[N-2]*x[1] + coeff[N-1]*x[2] + coeff[N]*x[3]
        for i = 1:N
            f += coeff[i]*rbf(x, centers[i])
        end
        return f
    end
    points, faces = isosurface(f_eval, MarchingTetrahedra(), origin=origin, widths=dims, samples=(pts,pts,pts))
    return points, faces
end

""" Construct mesh compatible with MeshCat
Args (output of isosurface):
    -points: Array of vertices
    -faces: Vertex connectivity graph
"""
function make_mesh(points, faces)
    npoints = size(points,1)
    nfaces = size(faces,1)
    meshpoints = Vector{Point3f0}(undef,npoints)
    for i = 1:npoints
        meshpoints[i] = Point3f0(points[i][1],points[i][2],points[i][3])
    end
    meshfaces = Vector{GLTriangleFace}(undef,nfaces)
    for i = 1:nfaces
        meshfaces[i] = TriangleFace(faces[i][1],faces[i][2],faces[i][3])
    end
    return Mesh(meshpoints, meshfaces)
end

""" Make implicit surface mesh with sdf """
function make_mesh_sdf(coeff::SVector{M,T}, centers::SVector{N, SVector{P,T}}, origin::SVector{3,T}, dims::SVector{3,T}, resolution::T) where {M,N,P,T}
    function f_eval(x)
        f = coeff[M-3]
        f += SVector{3,T}(coeff[M-2], coeff[M-1], coeff[M])'*x
        for i = 1:N
            f += coeff[i]*rbf(x, centers[i])
        end
        return f
    end
    sdf = gt.SignedDistanceField(f_eval, gt.HyperRectangle(gt.Vec(origin), gt.Vec(dims)), resolution)
    mesh = gt.HomogenousMesh(sdf, MarchingTetrahedra())
    return mesh
end
