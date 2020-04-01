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

---

weights = projection_weights(simplex)  # return deltas[:,i] ./ sum(deltas[:,i]) -> array
min_weight, index_to_replace = findmin(weights)
best_point = linear_combination(weights, simplex)	

---

since findmin(weights) is called, then set the elements of the weights to be 0 and then it will attempt to replace them.
