# """ Fill 1D boundaries: transimissive """
# function fill_trans_bc(g::Grid)
#     for k = 1:3, i = 1:g.n_g
#         g.cc_cons[i, k] = g.cc_cons[g.jlo, k]
#         g.cc_cons[g.jhi+i, k] = g.cc_cons[g.jhi, k]
#     end
# end


""" Fill 2D boundaries: transimissive """
function fill_trans_bc(g::Grid)
  # for i = 1:g.n_g
  #   g.cc_cons[i, :, :] .= g.cc_cons[g.i_x1_a, :, :]
  #   g.cc_cons[g.i_x1_b+i, :, :] .= g.cc_cons[g.i_x1_b, :, :]
  #   g.cc_cons[:, i, :] .= g.cc_cons[:, g.i_x2_a, :]
  #   g.cc_cons[:, g.i_x2_b+i, :] .= g.cc_cons[:, g.i_x2_b, :]
  # end
  @views @. g.cc_cons[g.i_x1g_a:g.i_x1g_b, :, :] = g.cc_cons[g.i_x1_a:g.i_x1_a, :, :]
  @views @. g.cc_cons[g.i_x1g_c:g.i_x1g_d, :, :] = g.cc_cons[g.i_x1_b:g.i_x1_b, :, :]
  @views @. g.cc_cons[:, g.i_x2g_a:g.i_x2g_b, :] = g.cc_cons[:, g.i_x2_a:g.i_x2_a, :]
  @views @. g.cc_cons[:, g.i_x2g_c:g.i_x2g_d, :] = g.cc_cons[:, g.i_x2_b:g.i_x2_b, :]
end


""" Fill periodic boundary conditions in 2D space """
function fill_periodic_bc(g::Grid)
  # for k = 1:4, i = 1:g.n_g
  #     g.cc_cons[i, :, k] = g.cc_cons[i_x1_b-g.n_g+i, :, k]
  #     g.cc_cons[i_x1_b+i, :, k] = g.cc_cons[g.i_x1_a+i-1, :, k]
  #     g.cc_cons[:, i, k] = g.cc_cons[:, g.i_x2_b-g.n_g+i, k]
  #     g.cc_cons[:, g.i_x2_b+i, k] = g.cc_cons[:, g.i_x2_a+i-1, k]
  # end
  for i = 1:g.n_g
    @. g.cc_cons[i, :, :] = g.cc_cons[i_x1_b-g.n_g+i, :, :]
    @. g.cc_cons[i_x1_b+i, :, :] = g.cc_cons[g.i_x1_a+i-1, :, :]
    @. g.cc_cons[:, i, :] = g.cc_cons[:, g.i_x2_b-g.n_g+i, :]
    @. g.cc_cons[:, g.i_x2_b+i, :] = g.cc_cons[:, g.i_x2_a+i-1, :]
  end
  # TODO: - finish vectorizing
  # @. g.cc_cons[g.i_x1g_a:g.i_x1g_b, :, :] = g.cc_cons[g.i_x1_b-g.n_g+1:i-g.i_x1_b, :, :]
  # @. g.cc_cons[i_x1_b+i, :, :] = g.cc_cons[g.i_x1_a+i-1, :, :]
  # @. g.cc_cons[:, i, :] = g.cc_cons[:, g.i_x2_b-g.n_g+i, :]
  # @. g.cc_cons[:, g.i_x2_b+i, :] = g.cc_cons[:, g.i_x2_a+i-1, :]
end
