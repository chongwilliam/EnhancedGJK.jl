""" Construct variational implicit surface given unordered set of points in R3 """

function rbf(x, y)
    return norm(x-y)^3
end

function A_fun!(A::MMatrix{M,N,T}, p::MMatrix{P,Q,T}) where {M,N,P,Q,T}
    """ Form interpolation matrix """
    for i = 1:P
        A[i,P+2:N] = p[i,:]
        A[P+2:M,i] = p[i,:]
        for j = 1:P
            A[i,j] = rbf(p[i,:],p[j,:])
        end
    end
    A[P+1,1:P] .= 1
    A[1:P,P+1] .= 1
end

function h_fun!(h::MVector{M,T}, f::MVector{N,T}) where {M,N,T}
    """ Form h matrix corresponding to p """
    h[1:N] = f[1:N]
end

function implicit_solve!(c::MVector{P,T}, A::MMatrix{M,N,T}, h::MVector{O,T}) where {P,M,N,O,T}
    """ Solve implicit surface interpolation equation with gjk distance constraints """
    c[:] = A \ h
end

function eval_surface(x::SVector{3,T}, c::MVector{M,T}, p::MMatrix{P,Q,T}) where {M,P,Q,T}
    """ Evaluate implicit surface value given point x (single evaluation) """
    f = c[M-3] + c[M-2]*x[1] + c[M-1]*x[2] + c[M]*x[3]
    for i = 1:M-4
        f += c[i]*rbf(x, p[i,:])
    end
    return f
end

# # Meshing functions
# function make_geom(c, center, dims, pts)
#     function f_eval(x)
#         f = c[n_pts+4-3] + c[n_pts+4-2]*x[1] + c[n_pts+4-1]*x[2] + c[n_pts+4]*x[3]
#         for i = 1:n_pts
#             f += c[i]*rbf(x,p[i,:])
#         end
#         return f
#     end
#     points, faces = isosurface(f_eval, MarchingCubes(), origin=center, widths=dims, samples=(pts,pts,pts))
#     return points, faces
# end
#
# function make_mesh(points, faces)
#     """ Construct mesh compatible with MeshCat
#     Args (output of isosurface):
#         -points: Array of vertices
#         -faces: Vertex connectivity graph
#     """
#     npoints = size(points,1)
#     nfaces = size(faces,1)
#     meshpoints = Vector{Point3f0}(undef,npoints)
#     for i = 1:npoints
#         meshpoints[i] = Point3f0(points[i][1],points[i][2],points[i][3])
#     end
#     meshfaces = Vector{GLTriangleFace}(undef,nfaces)
#     for i = 1:nfaces
#         meshfaces[i] = TriangleFace(faces[i][1],faces[i][2],faces[i][3])
#     end
#     return GeometryBasics.Mesh(meshpoints, meshfaces)
# end
