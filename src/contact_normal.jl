"""
Calculates the contact normal for each point that is in contact with the privileged object from GJK.
The collection of neighbor points for each contact point is determined from the active simplex points from
GJK results of the neighbor points matching with the point in contact.  Normal and tangential directions are computed
from PCA of the covariance of the set of points: cov = sum(y in neighbor set) outer((y-mean), (y-mean))
Arguments:
result: GJKResult of point to be investigated
points: list of GJKResult of points in point cloud
min_points: minimum number of points to be considered valid contact approximation
tol: distance tolerance for neighboring points
"""
function normal(result::GJKResult, points, min_points, tol, opt)  # neighbor points is a list of other GJKResult(s)
    if result.in_collision || separation_distance(result) > tol  # return null if in collision or too far away
        return -1
    else
        mean = MVector{3,Float64}(zeros(3,1))
        neighbor_points = []
        normal_vector = MVector{3,Float64}(zeros(3,1))  # averaged direction from mesh to point
        min_pts_flag = 0  # flag if minimum points is obtained

        # Add the current result as the first point
        if opt == 1  # option to use the closest point on the mesh from neighbor points
            push!(neighbor_points, result.poseA(result.closest_point_in_body.a))  # use projection onto mesh
            mean += result.poseA(result.closest_point_in_body.a)
            normal_vector += result.poseB(result.closest_point_in_body.b) - result.poseA(result.closest_point_in_body.a)
            num_points = 1
        else  # option to use the neighbor points directly
            push!(neighbor_points, result.poseB(result.closest_point_in_body.b))  # use environment point itself
            mean += result.poseB(result.closest_point_in_body.b)
            normal_vector += result.poseB(result.closest_point_in_body.b) - result.poseA(result.closest_point_in_body.a)
            num_points = 1
        end

        # Check if neighbor point's active simplex set is a subset of the main result's active simplex set
        for i = 1:length(points)
            if points[i].in_collision || separation_distance(points[i]) > tol  # ignore collision points or points that are too far away
                continue
            else
                cnt = 0
                M = length(result.mesh_simplex)
                N = length(points[i].mesh_simplex)
                # Neighbor point active simplex set subset check
                for j = 1:M  # main result
                    for k = 1:N  # potential neighbor result
                        if points[i].poseA(points[i].mesh_simplex[k]) == result.poseA(result.mesh_simplex[j])
                            cnt += 1
                            break
                        end
                    end
                end
                if (cnt == M) || (cnt == N && N < M) # enforce subset of simplex
                # if cnt == M  # enforce exactly matching simplex
                    if opt == 1
                        if points[i].poseA(points[i].closest_point_in_body.a) in neighbor_points  # check if already in the list of neighbor points
                            continue
                        else
                            push!(neighbor_points, points[i].poseA(points[i].closest_point_in_body.a))  # using projection onto face
                            mean += points[i].poseA(points[i].closest_point_in_body.a)  # using projection onto face
                            normal_vector += points[i].poseB(points[i].closest_point_in_body.b) - points[i].poseA(points[i].closest_point_in_body.a)
                            num_points += 1
                        end
                    else
                        if points[i].poseB(points[i].closest_point_in_body.b) in neighbor_points
                            continue
                        else
                            push!(neighbor_points, points[i].poseB(points[i].closest_point_in_body.b))  # using point itself
                            mean += points[i].poseB(points[i].closest_point_in_body.b)  # using point itself
                            normal_vector += points[i].poseB(points[i].closest_point_in_body.b) - points[i].poseA(points[i].closest_point_in_body.a)
                            num_points += 1
                        end
                    end
                end
                # Termination check
                if num_points >= min_points
                    # mean = mean / num_points
                    # normal_vector = normal_vector / num_points
                    min_pts_flag = 1
                    break  # break if want exact number of neighbor points = min points, otherwise comment out
                end
            end
        end

        # Calculate contact directions
        if min_pts_flag == 1
            mean = mean / num_points
            normal_vector = normal_vector / num_points
            cov = MMatrix{3,3,Float64}(zeros(3,3))
            for i = 1:num_points
                cov += (neighbor_points[i]-mean) * (neighbor_points[i]-mean)'  # covariance matrix
            end
            evalue, evector = eigen(cov)
            evector_sort = evector[:, sortperm(evalue)]  # sort from smallest to largest
            evector_final = []
            # Get normal vector
            vector = evector_sort[:,1]
            # Correct direction of normal to point towards the environment
            if dot(vector/LinearAlgebra.norm(vector), normal_vector/LinearAlgebra.norm(normal_vector)) < 0
                push!(evector_final, -vector)
            else
                push!(evector_final, vector)
            end
            push!(evector_final, evector_sort[:,2])
            push!(evector_final, evector_sort[:,3])
            return [mean, evector_final, num_points, neighbor_points]  # mean point, contact directions, num_points, list of neighbor points
        else
            return -1  # not enough points flag
        end
    end
end

"""
Unit test for the general algorithm.
"""
function normal_test(points)
    mean = MVector{3,Float64}(zeros(3,1))
    for i = 1:length(points)
        mean += points[i]
    end
    mean = mean / length(points)
    cov = MMatrix{3,3,Float64}(zeros(3,3))
    for i = 1:length(points)
        cov += (points[i]-mean) * (points[i]-mean)'
    end
    evalue, evector = eigen(cov)  # eigenvalues/vectors not ordered
    return evector
end
