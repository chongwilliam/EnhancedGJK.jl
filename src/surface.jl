""" Construct variational implicit surface given unordered set of points in R3 """

""" Duchon's interpolant """
function rbf(x::N, y::T) where {N,T}
    return norm(x-y)^3
end

""" Form interpolation matrix """
function A_fun!(A::MMatrix{P,P,T}, centers::SVector{M, SVector{N,T}}) where {M,N,P,T}
    for i = 1:M
        A[i,P-2:P] = centers[i]
        A[P-2:P,i] = centers[i]
        for j = 1:M
            A[i,j] = rbf(centers[i],centers[j])
        end
    end
    A[M+1,1:M] .= 1
    A[1:M,M+1] .= 1
    A[M+1:P,M+1:P] .= 0
end

""" Form b vector corresponding to signed distance values """
function b_fun!(b::MVector{M,T}, distances::SVector{N,T}) where {M,N,T}
    for i = 1:N
        b[i] = distances[i]
    end
    b[N+1:N+4] .= 0
end

""" Evaluate implicit surface value given point x (single evaluation) """
function eval_surface(c::SVector{M,T}, centers::SVector{N, SVector{P,T}}, x::SVector{P,T}) where {M,N,P,T}
    # f = c[M-3] + c[M-2]*x[1] + c[M-1]*x[2] + c[M]*x[3]
    f = c[M-3]
    f += SVector{3,T}(c[M-2], c[M-1], c[M])'*x
    for i = 1:N
        f += c[i]*rbf(x, centers[i])
    end
    return f
end

# """ Evalute implicit surface value given point x for ForwardDiff. x = [x; coeff; centers] """
# function eval_surface(x::Vector)
#     N = length(x)  # (3; M+4; M) = 2M + 7
#     M = Int((N-7)/2)
#     z = x[1:3]
#     c = x[4:M+4]
#     centers = x[M+5:end]
#     f = c[M+4]
#     f += SVector{3,T}(c[M+3], c[M+2], c[M+1])'*z
#     for i = 1:M
#         f += c[i]*rbf(z, centers[i])
#     end
#     return f
# end
#
# function distance_gradient(c::SVector{M,T}, centers::SVector{N, SVector{P,T}}, x::SVector{P,T}) where {M,N,P,T}
#     input = Vector{T}(undef, 3+M+N)
#     input[1:3] = x
#     input[4:M+4] = c
#     input[M+5:end] = centers
#     ∇f = ForwardDiff.gradient(eval_surface, input)
#     return ∇f
# end


""" Construct implicit surface from closest N environment points to a face """
function construct_surface(centers::SVector{M, SVector{N,T}}, distances::SVector{M,T}) where {M,N,T}
    A = MMatrix{M+4, M+4, T}(undef)
    b = MVector{M+4, T}(undef)
    A_fun!(A, centers)
    b_fun!(b, distances)
    c = A \ b  # get coefficients
    return SVector{M+4,Float64}(c)
end

""" Create vector for rbf centers """
function construct_centers(verts::Vector{GeometryTypes.Point{3,T}}, points::SVector{M,SVector{3,T}}, distances::SVector{M,T}) where {M,T,U}
    N = length(verts)
    centers = MVector{M+N, SVector{3,Float64}}(undef)
    dist = MVector{M+N, Float64}(zeros(M+N))
    for i = 1:N
        centers[i] = SVector{3,Float64}(verts[i])
    end
    for i = 1:M
        centers[N+i] = points[i]
        dist[N+i] = distances[i]
    end
    return SVector(centers), SVector(dist)
end

# Meshing functions
""" Return points and faces of the isosurface of the implicit surface """
function make_geom(c::SVector{M,T}, centers::SVector{N, SVector{P,T}}, origin::SVector{3,T}, dims::SVector{3,T}, pts::Int) where {M,N,P,T}
    function f_eval(x)
        f = c[N-3] + c[N-2]*x[1] + c[N-1]*x[2] + c[N]*x[3]
        for i = 1:N
            f += c[i]*rbf(x,centers[i])
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
function make_mesh_sdf(c::SVector{M,T}, centers::SVector{N, SVector{P,T}}, origin::SVector{3,T}, dims::SVector{3,T}, resolution::T) where {M,N,P,T}
    function f_eval(x)
        f = c[M-3]
        f += SVector{3,T}(c[M-2], c[M-1], c[M])'*x
        for i = 1:N
            f += c[i]*rbf(x, centers[i])
        end
        return f
    end
    sdf = gt.SignedDistanceField(f_eval, gt.HyperRectangle(gt.Vec(origin), gt.Vec(dims)), resolution)
    mesh = gt.HomogenousMesh(sdf, MarchingTetrahedra())
    return mesh
end
