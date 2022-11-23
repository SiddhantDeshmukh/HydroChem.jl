# Requires grid.jl, flux.jl

function minmod(x, y, z)
  0.25 * abs(sign(x) + sign(y)) * (sign(x) + sign(z)) *
  min(abs(x), abs(y), abs(z))
end



function reconstruct1st(g::Grid, axis::Int)
  if axis == 1
    for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
      g.prim_L[i, j, k] = g.cc_prim[i, j, k]
      g.prim_R[i, j, k] = g.cc_prim[i+1, j, k]
    end
    for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
      g.cs_L[i, j] = g.cs[i, j]
      g.cs_R[i, j] = g.cs[i+1, j]
    end
  elseif axis == 2
    for k = 1:g.nvars, j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
      g.prim_L[i, j, k] = g.cc_prim[i, j, k]
      g.prim_R[i, j, k] = g.cc_prim[i, j+1, k]
    end
    for j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
      g.cs_L[i, j] = g.cs[i, j]
      g.cs_R[i, j] = g.cs[i, j+1]
    end
  end
end


""" Interpolate primitive variables (ρ, vx, vy, p, cs) in the x and y components """
function reconstruct2nd(g::Grid, axis::Int; theta::Float64=1.5)
  @debug "Inside reconstruct2nd"
  cdiff = zero(g.cc_prim)
  if axis == 1
    for k = 1:g.nvars  # reconstruct everything, including passive scalars
      c = @view g.cc_prim[:, :, k]     # c = rho, vx, vy, vz, pressure for k = 1:5
      # @show k, c[g.mid_x1, g.mid_x2]
      for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b+1
        cdiff[i, j, k] = 0.5 * minmod(theta * (c[i, j] - c[i-1, j]),
          0.5 * (c[i+1, j] - c[i-1, j]),
          theta * (c[i+1, j] - c[i, j]))
      end
    end
    for i = g.i_x1_a-1:g.i_x1_b
      @. g.prim_L[i, g.i_x2_a:g.i_x2_b, :] = g.cc_prim[i, g.i_x2_a:g.i_x2_b, :] + cdiff[i, g.i_x2_a:g.i_x2_b, :]
      @. g.prim_R[i, g.i_x2_a:g.i_x2_b, :] = g.cc_prim[i+1, g.i_x2_a:g.i_x2_b, :] - cdiff[i+1, g.i_x2_a:g.i_x2_b, :]
    end
    # @show g.rho_L[g.mid_x1, g.mid_x2, :], g.p_L[g.mid_x1, g.mid_x2, :], g.rho_R[g.mid_x1, g.mid_x2, :], g.p_R[g.mid_x1, g.mid_x2, :]
    # @show g.prim_L[g.mid_x1, g.mid_x2, :]
    # @show g.prim_R[g.mid_x1, g.mid_x2, :]
    # @show cdiff[g.mid_x1, g.mid_x2, :]
    # @show g.cc_prim[g.mid_x1, g.mid_x2, :]
    # exit()
    for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
      g.cs_L[i, j] = prim2cs(g.rho_L[i, j], g.p_L[i, j], g.gamma)
      g.cs_R[i, j] = prim2cs(g.rho_R[i, j], g.p_R[i, j], g.gamma)
    end
  elseif axis == 2
    for k = 1:g.nvars
      c = @view g.cc_prim[:, :, k]     # c = rho, vx, vy, vz, pressure for k = 1:5
      for j = g.i_x2_a-1:g.i_x2_b+1, i = g.i_x1_a:g.i_x1_b
        cdiff[i, j, k] = 0.5 * minmod(theta * (c[i, j] - c[i, j-1]),
          0.5 * (c[i, j+1] - c[i, j-1]),
          theta * (c[i, j+1] - c[i, j]))
      end
    end
    for k = 1:g.nvars, j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
      g.prim_L[i, j, k] = g.cc_prim[i, j, k] + cdiff[i, j, k]
      g.prim_R[i, j, k] = g.cc_prim[i, j+1, k] - cdiff[i, j+1, k]
    end
    for j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
      g.cs_L[i, j] = prim2cs(g.rho_L[i, j], g.p_L[i, j], g.gamma)
      g.cs_R[i, j] = prim2cs(g.rho_R[i, j], g.p_R[i, j], g.gamma)
    end
  end
end