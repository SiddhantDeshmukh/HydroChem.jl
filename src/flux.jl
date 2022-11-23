""" 2d: From primitive variables (rho, vel, pressure) to conservative variables """
function prim2cons!(ρ, v_x1, v_x2, v_x3, p, cons::Array{Float64,3}, gamma::Float64)
  @. cons[:, :, 1] = ρ  # density -> density
  @. cons[:, :, 2] = ρ * v_x1  # velocity -> mass flux
  @. cons[:, :, 3] = ρ * v_x2  # velocity -> mass flux
  @. cons[:, :, 4] = ρ * v_x3  # velocity -> mass flux
  @. cons[:, :, 5] = p / (gamma - 1) + 0.5 * ρ * (v_x1^2 + v_x2^2 + v_x3^2)  # pressure -> energy
  return cons
end

function prim2cons!(g::Grid; convert_centres=true, convert_interfaces=false)
  if convert_centres
    prim2cons!(g.rho, g.v_x1, g.v_x2, g.v_x3, g.p, g.cc_cons, g.gamma)
  end
  if convert_interfaces
    prim2cons!(g.rho_L, g.v_x1_L, g.v_x2_L, g.v_x3_L, g.p_L, g.cons_L, g.gamma)
    prim2cons!(g.rho_R, g.v_x1_R, g.v_x2_R, g.v_x3_L, g.p_R, g.cons_R, g.gamma)
  end
  return g
end


""" From prims to sound speed """
function prim2cs(ρ::Float64, p::Float64, γ::Float64)
  if (ρ < 0) || (p < 0)
    @show γ, ρ, p
  end

  return sqrt(γ * p / ρ)
end

""" 2d: Reconstruct rho, vel, E, epsilon, and pressure from g.cc_cons """
function cons2prim!(g::Grid)
  @. g.rho = g.cc_cons[:, :, 1]
  @. g.v_x1 = g.cc_cons[:, :, 2] / g.rho
  @. g.v_x2 = g.cc_cons[:, :, 3] / g.rho
  @. g.v_x3 = g.cc_cons[:, :, 4] / g.rho
  @. g.E = g.cc_cons[:, :, 5]
  @. g.epsilon = g.E / g.rho - 0.5 * (g.v_x1^2 + g.v_x2^2)
  @. g.p = g.lambda * g.rho * g.epsilon
  @. g.cs = sqrt(g.gamma * g.p / g.rho)

  return g
end


""" 2D: Calculate flux from prims (rho, vel, pressure) """
function calc_flux!(ρ, v_x1, v_x2, v_x3, p, γ::Float64, f_cons_x1, f_cons_x2;
  ps=nothing)
  @. begin
    # f_cons_x1[:, :, g.i_rho] = ρ * v_x1
    # f_cons_x1[:, :, g.i_rhoun] = ρ * v_x1^2 + p
    # f_cons_x1[:, :, g.i_rhout] = ρ * v_x1 * v_x2
    # f_cons_x1[:, :, g.i_rhoux3] = ρ * v_x1 * v_x3
    # f_cons_x1[:, :, g.i_rhouE] = (p / (γ - 1) + 0.5 * ρ * (v_x1^2 + v_x2^2) + p) * v_x1
    # f_cons_x2[:, :, g.i_rho] = ρ * v_x2
    # f_cons_x2[:, :, g.i_rhoun] = ρ * v_x1 * v_x2
    # f_cons_x2[:, :, g.i_rhout] = ρ * v_x2^2 + p
    # f_cons_x2[:, :, g.i_rhoux3] = ρ * v_x1 * v_x3
    # f_cons_x2[:, :, g.i_rhouE] = (p / (γ - 1) + 0.5 * ρ * (v_x1^2 + v_x2^2) + p) * v_x2
    f_cons_x1[:, :, 1] = ρ * v_x1
    f_cons_x1[:, :, 2] = ρ * v_x1^2 + p
    f_cons_x1[:, :, 3] = ρ * v_x1 * v_x2
    f_cons_x1[:, :, 4] = ρ * v_x1 * v_x3
    f_cons_x1[:, :, 5] = (p / (γ - 1) + 0.5 * ρ * (v_x1^2 + v_x2^2) + p) * v_x1
  end
  if !isnothing(ps)
    for k = axes(ps, 3)
      @. f_cons_x1[:, :, 5+k] = ps[:, :, k] * v_x1
    end
  end

  @. begin
    f_cons_x2[:, :, 1] = ρ * v_x2
    f_cons_x2[:, :, 2] = ρ * v_x1 * v_x2
    f_cons_x2[:, :, 3] = ρ * v_x2^2 + p
    f_cons_x2[:, :, 4] = ρ * v_x1 * v_x3
    f_cons_x2[:, :, 5] = (p / (γ - 1) + 0.5 * ρ * (v_x1^2 + v_x2^2) + p) * v_x2
  end
  if !isnothing(ps)
    for k = axes(ps, 3)
      @. f_cons_x2[:, :, 5+k] = ps[:, :, k] * v_x2
    end
  end
end

# function calc_flux!(ρ, v_x1, v_x2, v_x3, p, γ::Float64)
#   # f_cons_x1 = Array{Float64}(undef, size(ρ, 1), 4)
#   f_cons_x1 = Array{Float64}(undef, size(ρ, 1), size(ρ, 2), 4)
#   f_cons_x2 = Array{Float64}(undef, size(ρ, 1), size(ρ, 2), 4)
#   calc_flux!(ρ, v_x1, v_x2, v_x3, p, γ, f_cons_x1, f_cons_x2)
#   return f_cons_x1, f_cons_x2
# end

# function calc_flux!(g::Grid)
#   calc_flux!(g.rho, g.v_x1, g.v_x2, g.v_x3, g.p, g.gamma, g.f_cons_x1, g.f_cons_x2; ps=g.ps)
#   return g.f_cons_x1, g.f_cons_x2
# end
