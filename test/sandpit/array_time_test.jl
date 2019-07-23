# function timing_test(ho::Int64 = 40)
#     len = (ho+1) * (2ho+1)
#     # MM = Matrix{Float64}(undef,len,len)
#     MM = Matrix{Array{Int64}}(undef,len,len);
#     function fillMM!(MM)
#         ind1 = 1
#         # for l = 0:ho for m = -l:l, s1 = 1:S
#         for ds = 0:ho for dm = -ho:ho
#             ind2 = 1
#             # for dl = 0:ho for dm = -dl:dl, s2 = 1:S
#             for s = 0:ho for m = -ho:ho
#                 MM[ind1,ind2] = [s, m, ds, dm]
#                 ind2 += 1
#             end
#             end
#             ind1 += 1
#         end
#         end
#         return MM
#     end
#
#     function fillMM2!(MM)
#         # ind2 = 1
#         # for l = 0:ho for m = -l:l, s1 = 1:S
#         # for s = 0:ho, m = -ho:ho
#         #     # ind1 = 1
#         #     # for dl = 0:ho for dm = -dl:dl, s2 = 1:S
#         #     # MM[:,ind2] = [ [s, m, ds, dm] for dm = -ho:ho, ds = 0:ho][:]
#         #     MM[:,ind2] = vcat( [[[s, m, ds, dm] for dm = -ho:ho] for ds = 0:ho]... )
#         #     ind2 += 1
#         # end
#         MM[:] = hcat([
#             vcat( [[[s, m, ds, dm] for dm = -ho:ho] for ds = 0:ho]... )
#         for m = -ho:ho, s = 0:ho]... )
#         return MM
#     end
#     println("fillMM2")
#     @time fillMM2!(MM);
#     println("fillMM")
#     @time fillMM!(MM);
#     return nothing
# end

# reshape(MM[1,:], (2ho+1,ho+1))
#
# @time MM2 = reshape(
#         [ [s, m, ds, dm] for dm in -ho:ho, ds in 0:ho, m in -ho:ho, s in 0:ho]
# , ((2ho+1)*(ho+1), (2ho+1)*(ho+1)));
#
# using LinearAlgebra
#
# norm(norm.(MM2 - MM))
#
# # @time MM2 = reshape(
# #     [[s, m, ds, dm] for m in -ho:ho, s in 0:ho, dm in -ho:ho, ds in 0:ho]
# # , ((2ho+1)*(ho+1), (2ho+1)*(ho+1)))
#
# reshape(MM2[1,:], (2ho+1,ho+1))
#
# MM2 = [ [ds, dm] for dm in -ho:ho, ds in 0:ho]
#
# MM = deepcopy(MM2)
# for ds = 0:ho for dm = -ho:ho
#     MM[dm+ho+1,ds+1] = [ds, dm]
# end end
#
# MM
