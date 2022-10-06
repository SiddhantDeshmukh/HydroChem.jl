# Requires grid.jl, flux.jl, reconstruction.jl
# Riemann Solvers

# uL - aL, uR - aR
function util_speed_davis1(uL, uR, csL, csR)
  return uL - csL, uR - csR
end


# min(uL - aL, uR - aR), max(uL + aL, uR + aR)
function util_speed_davis2(uL, uR, csL, csR)
  return min(uL - csL, uR - csR), max(uL + csL, uR + csR)
end


# Wave speed estimater: Roe 1981 JCP 43:357–372.
# need L and R states of rho, v, E, p
function util_speed_roe(rho_L, vL, EL, pL, rho_R, vR, ER, pR, gamma)
  uhat = (sqrt(rho_L) * vL + sqrt(rho_R) * vR) / (sqrt(rho_L) + sqrt(rho_R))
  HL = (EL + pL) / rho_L
  HR = (ER + pR) / rho_R
  Hhat = (sqrt(rho_L) * HL + sqrt(rho_R) * HR) / (sqrt(rho_L) + sqrt(rho_R))
  if (Hhat < 0) || (uhat < 0)
    @show Hhat, uhat
  end
  # @show Hhat, uhat
  if (Hhat < 0.5 * uhat^2)
    @show Hhat, uhat
    @show rho_L, vL, EL, pL
    @show rho_R, vR, ER, pR
  end
  # Following line gives sqrt of negative for sod shock test
  # i.e. 0.5 * uhat^2 > Hhat for a point
  ahat = sqrt((gamma - 1) * (Hhat - 0.5 * uhat^2))
  return uhat - ahat, uhat + ahat
end



""" Second-order HLL. This should not change g.cc_cons """
function hll(g::Grid, reconstruct::Function=reconstruct2nd)
  println("Inside HLL solver")
  # x component
  # interpolate_x(g)
  reconstruct(g, 1)
  println("Done reconstruction on x1")
  # no potential sqrt error
  prim2cons!(g.rho_L, g.v_x1_L, g.v_x2_L, g.v_x3_L, g.p_L, g.cons_L, g.gamma)
  prim2cons!(g.rho_R, g.v_x1_R, g.v_x2_R, g.v_x3_R, g.p_R, g.cons_R, g.gamma)
  # alpha_plus = max.(0.0, g.v_x1_L .+ g.csL, g.v_x1_R .+ g.csR)
  # alpha_minus = max.(0.0, -g.v_x1_L .+ g.csL, -g.v_x1_R .+ g.csR)
  lam_minus = zero(g.cs)
  lam_plus = zero(g.cs)
  for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
    # # davis2
    # lam_minus[i, j], lam_plus[i, j] = util_speed_davis2(
    #     g.v_x1_L[i, j], g.v_x1_R[i, j], g.csL[i, j], g.csR[i, j])
    # roe
    lam_minus[i, j], lam_plus[i, j] = util_speed_roe(
      g.rho_L[i, j], g.v_x1_L[i, j], g.cons_L[i, j, 5], g.p_L[i, j],
      g.rho_R[i, j], g.v_x1_R[i, j], g.cons_R[i, j, 5], g.p_R[i, j],
      g.gamma)
  end
  alpha_minus = max.(0.0, -1.0 .* lam_minus)
  alpha_plus = max.(0.0, lam_plus)
  calc_flux!(g.rho_L, g.v_x1_L, g.v_x2_L, g.v_x3_L, g.p_L, g.gamma, g.f_cons_x1, g.f_cons_x2)
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
    g.fhll[i, j, k] = alpha_plus[i, j] * g.f_cons_x1[i, j, k]
  end
  calc_flux!(g.rho_R, g.v_x1_R, g.v_x2_R, g.v_x3_R, g.p_R, g.gamma, g.f_cons_x1, g.f_cons_x2)
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
    g.fhll[i, j, k] = (g.fhll[i, j, k] + alpha_minus[i, j] * g.f_cons_x1[i, j, k]
                       -
                       alpha_plus[i, j] * alpha_minus[i, j] *
                       (g.cons_R[i, j, k] - g.cons_L[i, j, k])) /
                      (alpha_plus[i, j] + alpha_minus[i, j])
  end
  # calculate L(u) contributed by the x component of the flux
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    g.lu[i, j, k] = -(g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx1
  end

  # y component
  # interpolate_y(g)
  reconstruct(g, 2)
  println("Done reconstruction on x2")
  # no potential sqrt error
  prim2cons!(g.rho_L, g.v_x1_L, g.v_x2_L, g.v_x3_L, g.p_L, g.cons_L, g.gamma)
  prim2cons!(g.rho_R, g.v_x1_R, g.v_x2_R, g.v_x3_R, g.p_R, g.cons_R, g.gamma)
  # alpha_plus = max.(0.0, g.v_x1_L .+ g.csL, g.v_x1_R .+ g.csR)
  # alpha_minus = max.(0.0, -g.v_x1_L .+ g.csL, -g.v_x1_R .+ g.csR)
  lam_minus = zero(g.cs)
  lam_plus = zero(g.cs)
  for j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    # # davis2
    # lam_minus[i, j], lam_plus[i, j] = util_speed_davis2(
    #     g.v_x2_L[i, j], g.v_x2_R[i, j], g.csL[i, j], g.csR[i, j])
    # roe
    lam_minus[i, j], lam_plus[i, j] = util_speed_roe(
      g.rho_L[i, j], g.v_x2_L[i, j], g.cons_L[i, j, 5], g.p_L[i, j],
      g.rho_R[i, j], g.v_x2_R[i, j], g.cons_R[i, j, 5], g.p_R[i, j],
      g.gamma)
  end
  alpha_minus = max.(0.0, -1.0 .* lam_minus)
  alpha_plus = max.(0.0, lam_plus)

  calc_flux!(g.rho_L, g.v_x1_L, g.v_x2_L, g.v_x3_L, g.p_L, g.gamma, g.f_cons_x1, g.f_cons_x2)
  for k = 1:g.nvars, j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    g.fhll[i, j, k] = alpha_plus[i, j] * g.f_cons_x2[i, j, k]
  end
  calc_flux!(g.rho_R, g.v_x1_R, g.v_x2_R, g.v_x3_R, g.p_R, g.gamma, g.f_cons_x1, g.f_cons_x2)
  for k = 1:g.nvars, j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    g.fhll[i, j, k] = (g.fhll[i, j, k] +
                       alpha_minus[i, j] * g.f_cons_x2[i, j, k] -
                       alpha_plus[i, j] * alpha_minus[i, j] *
                       (g.cons_R[i, j, k] - g.cons_L[i, j, k])) /
                      (alpha_plus[i, j] + alpha_minus[i, j])
    # @. g.fhll[i, j, :] = (g.fhll[i, j, :] +
    #                    alpha_minus[i, j] * g.f_cons_x2[i, j, :] -
    #                    alpha_plus[i, j] * alpha_minus[i, j] *
    #                    (g.cons_R[i, j, :] - g.cons_L[i, j, :])) /
    #                   (alpha_plus[i, j] + alpha_minus[i, j])
  end
  # calculate L(u) contributed by the y component of the flux
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    g.lu[i, j, k] += -(g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dx2
  end
  return g.lu
end


""" HLLC solver. Work in progress. Some bugs remains to be fixed. """
function hllc_var1(g::Grid, reconstruct::Function=reconstruct2nd)
  # x component
  # interpolate_x(g)
  reconstruct(g, 1)
  # no potential sqrt error
  prim2cons!(g.rho_L, g.v_x1_L, g.v_x2_L, g.v_x3_L, g.p_L, g.cons_L, g.gamma)
  prim2cons!(g.rho_R, g.v_x1_R, g.v_x2_R, g.v_x3_L, g.p_R, g.cons_R, g.gamma)
  sl = similar(g.cs)
  sr = similar(g.cs)
  ss = similar(g.cs)
  for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
    # # davis2
    # sl[i, j], sr[i, j] = util_speed_davis2(
    #     g.v_x1_L[i, j], g.v_x1_R[i, j], g.csL[i, j], g.csR[i, j])
    # roe
    sl[i, j], sr[i, j] = util_speed_roe(
      g.rho_L[i, j], g.v_x1_L[i, j], g.cons_L[i, j, 4], g.p_L[i, j],
      g.rho_R[i, j], g.v_x1_R[i, j], g.cons_R[i, j, 4], g.p_R[i, j],
      g.gamma)
    ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rho_L[i, j] * g.v_x1_L[i, j] *
                                          (sl[i, j] - g.v_x1_L[i, j]) - g.rho_R[i, j] * g.v_x1_R[i, j] *
                                                                        (sr[i, j] - g.v_x1_R[i, j])) / (g.rho_L[i, j] *
                                                                                                        (sl[i, j] - g.v_x1_L[i, j]) - g.rho_R[i, j] * (sr[i, j] - g.v_x1_R[i, j]))
  end
  calc_flux!(g.rho_L, g.v_x1_L, g.v_x2_L, g.p_L, g.gamma, g.f_cons_x1, g.f_cons_x2)
  fuL = copy(g.f_cons_x1)
  calc_flux!(g.rho_R, g.v_x1_R, g.v_x2_R, g.p_R, g.gamma, g.f_cons_x1, g.f_cons_x2)
  fuR = g.f_cons_x1
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a-1:g.i_x1_b
    fsl = ss[i, j] * (sl[i, j] * g.cons_L[i, j, k] - fuL[i, j, k])
    fsr = ss[i, j] * (sr[i, j] * g.cons_R[i, j, k] - fuR[i, j, k])
    if k == 2
      fsl += sl[i, j] * (g.pL[i, j] + g.rho_L[i, j] *
                                      (sl[i, j] - g.v_x1_L[i, j]) * (ss[i, j] - g.v_x1_L[i, j]))
      fsr += sr[i, j] * (g.pR[i, j] + g.rho_L[i, j] *
                                      (sr[i, j] - g.v_x1_R[i, j]) * (ss[i, j] - g.v_x1_R[i, j]))
    elseif k == 4
      fsl += sl[i, j] * (g.pL[i, j] + g.rho_L[i, j] *
                                      (sl[i, j] - g.v_x1_L[i, j]) * (ss[i, j] - g.v_x1_L[i, j])) * ss[i, j]
      fsr += sr[i, j] * (g.pR[i, j] + g.rho_L[i, j] *
                                      (sr[i, j] - g.v_x1_R[i, j]) * (ss[i, j] - g.v_x1_R[i, j])) * ss[i, j]
    end
    fsl /= sl[i, j] - ss[i, j]
    fsr /= sl[i, j] - ss[i, j]
    if 0 <= sl[i, j]
      g.fhll[i, j, k] = fuL[i, j, k]
    elseif sl[i, j] < 0 <= ss[i, j]
      g.fhll[i, j, k] = fsl
    elseif ss[i, j] < 0 <= sr[i, j]
      g.fhll[i, j, k] = fsr
    else
      g.fhll[i, j, k] = fuR[i, j, k]
    end
  end
  # calculate L(u)
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    g.lu[i, j, k] = -(g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx1
  end

  # y component
  # interpolate_y(g)
  reconstruct(g, 2)
  # no potential sqrt error
  prim2cons!(g.rho_L, g.v_x1_L, g.v_x2_L, g.p_L, g.cons_L, g.gamma)
  prim2cons!(g.rho_R, g.v_x1_R, g.v_x2_R, g.p_R, g.cons_R, g.gamma)
  sl = similar(g.cs)
  sr = similar(g.cs)
  ss = similar(g.cs)
  for j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    # # davis2
    # sl[i, j], sr[i, j] = util_speed_davis2(
    #     g.v_x2_L[i, j], g.v_x2_R[i, j], g.csL[i, j], g.csR[i, j])
    # roe
    sl[i, j], sr[i, j] = util_speed_roe(
      g.rho_L[i, j], g.v_x2_L[i, j], g.cons_L[i, j, 4], g.p_L[i, j],
      g.rho_R[i, j], g.v_x2_R[i, j], g.cons_R[i, j, 4], g.p_R[i, j],
      g.gamma)
    ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rho_L[i, j] * g.v_x2_L[i, j] *
                                          (sl[i, j] - g.v_x2_L[i, j]) - g.rho_R[i, j] * g.v_x2_R[i, j] *
                                                                        (sr[i, j] - g.v_x2_R[i, j])) / (g.rho_L[i, j] *
                                                                                                        (sl[i, j] - g.v_x2_L[i, j]) - g.rho_R[i, j] * (sr[i, j] - g.v_x2_R[i, j]))
  end
  calc_flux!(g.rho_L, g.v_x1_L, g.v_x2_L, g.p_L, g.gamma, g.f_cons_x1, g.f_cons_x2)
  fuL = copy(g.f_cons_x2)
  calc_flux!(g.rho_R, g.v_x1_R, g.v_x2_R, g.p_R, g.gamma, g.f_cons_x1, g.f_cons_x2)
  fuR = g.f_cons_x2
  for k = 1:g.nvars, j = g.i_x2_a-1:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    fsl = ss[i, j] * (sl[i, j] * g.cons_L[i, j, k] - fuL[i, j, k])
    fsr = ss[i, j] * (sr[i, j] * g.cons_R[i, j, k] - fuR[i, j, k])
    if k == 2
      fsl += sl[i, j] * (g.pL[i, j] + g.rho_L[i, j] *
                                      (sl[i, j] - g.v_x2_L[i, j]) * (ss[i, j] - g.v_x2_L[i, j]))
      fsr += sr[i, j] * (g.pR[i, j] + g.rho_L[i, j] *
                                      (sr[i, j] - g.v_x2_R[i, j]) * (ss[i, j] - g.v_x2_R[i, j]))
    elseif k == 4
      fsl += sl[i, j] * (g.pL[i, j] + g.rho_L[i, j] *
                                      (sl[i, j] - g.v_x2_L[i, j]) * (ss[i, j] - g.v_x2_L[i, j])) * ss[i, j]
      fsr += sr[i, j] * (g.pR[i, j] + g.rho_L[i, j] *
                                      (sr[i, j] - g.v_x2_R[i, j]) * (ss[i, j] - g.v_x2_R[i, j])) * ss[i, j]
    end
    fsl /= sl[i, j] - ss[i, j]
    fsr /= sl[i, j] - ss[i, j]
    if 0 <= sl[i, j]
      g.fhll[i, j, k] = fuL[i, j, k]
    elseif sl[i, j] < 0 <= ss[i, j]
      g.fhll[i, j, k] = fsl
    elseif ss[i, j] < 0 <= sr[i, j]
      g.fhll[i, j, k] = fsr
    else
      g.fhll[i, j, k] = fuR[i, j, k]
    end
  end
  # calculate L(u) contributed by the y component of the flux
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    g.lu[i, j, k] += -(g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dx2
  end
  return g.lu
end


# LAX scheme, 2D
function lax(g::Grid)
  F, G = calc_flux(g)
  for k = 1:g.nvars, j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
    # calculate L(u)
    g.lu[i, j, k] = -0.5 / g.dx1 * (F[i+1, j, k] - F[i-1, j, k]) -
                    0.5 / g.dx2 * (G[i, j+1, k] - G[i, j-1, k])
    # replace g.cc_cons with (g.cc_cons[i-1, j] + g.cc_cons[i+1, j]) / 2
    g.cc_cons[i, j, k] = 0.25 * (g.cc_cons[i-1, j, k] + g.cc_cons[i+1, j, k] +
                                 g.cc_cons[i, j-1, k] + g.cc_cons[i, j+1, k])
  end
  return g.lu
end
