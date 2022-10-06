using Parameters

@with_kw mutable struct Grid
  # 2D finite-volume Grid to hold all quantities
  n_x1::Int = 10
  n_x2::Int = 1
  n_g::Int = 2
  x1_len::Int = n_x1 + 2 * n_g
  x2_len::Int = n_x2 + 2 * n_g
  mid_x1::Int64 = floor(Int64, x1_len / 2)
  mid_x2::Int64 = floor(Int64, x2_len / 2)
  n_ps::Int64 = 0  # number of passive scalars
  nvars::Int = 5 + n_ps # (rho, ux1, ux2, ux3, p) + passive scalars

  # Spatial domain (without ghost cells)
  min_x1::Float64 = 0.0
  max_x1::Float64 = 1.0
  min_x2::Float64 = 0.0
  max_x2::Float64 = 1.0
  t::Float64 = 0.0
  gamma::Float64 = 1.4
  lambda::Float64 = gamma - 1

  # Grid spacing
  dx1::Float64 = (max_x1 - min_x1) / n_x1
  dx2::Float64 = (max_x2 - min_x2) / n_x2

  # Cell edges and centres
  x1_l::Vector{Float64} = (collect(0:x1_len-1) .- n_g) .* dx1 # left edges
  x1::Vector{Float64} = x1_l .+ dx1 / 2                       # centers
  x2_l::Vector{Float64} = (collect(0:x2_len-1) .- n_g) .* dx2 # left edges
  x2::Vector{Float64} = x2_l .+ dx2 / 2                       # centers

  # Cell-centred quantities
  # Conserved: (rho, rhoux1, rhoux2, rhoux3, rhoE)
  cc_cons::Array{Float64,3} = zeros(Float64, x1_len, x2_len, nvars)
  # Primitive: (rho, ux1, ux2, ux3, p)
  cc_prim::Array{Float64,3} = zeros(Float64, x1_len, x2_len, nvars)
  # Sound speed, energy, epsilon
  cs = zeros(x1_len, x2_len)
  E = zero(cs)
  epsilon = zero(cs)

  # Left and right interface
  cons_L::Array{Float64,3} = zero(cc_cons)
  cons_R::Array{Float64,3} = zero(cc_cons)
  prim_L::Array{Float64,3} = zero(cc_prim)
  prim_R::Array{Float64,3} = zero(cc_prim)
  cs_L = zero(cs)
  cs_R = zero(cs)

  # Fluxes
  f_cons_x1::Array{Float64,3} = zero(cc_cons)
  f_cons_x2::Array{Float64,3} = zero(cc_cons)

  # xfu = similar(u)
  # yfu = similar(u)
  lu = zero(cc_cons)  # what is this?
  fhll = zero(cc_cons)

  # Indices (to prevent ambiguity)
  # Spatial domain: the specifiers 'a', 'b' define the first and
  # last indices (on a Cartesian grid, these are [left, right] for x and
  # [lower, upper] for y)
  # For ghost cells, 'c' and 'd' are also defined such that the lower ghost
  # cells occupy ['a', 'b'] and the upper ones occupy ['c', 'd']
  # TODO:
  # - diagram describing this
  # x1 inside domain
  i_x1_a::Int = 1 + n_g
  i_x1_b::Int = n_x1 + n_g
  # x1 ghost cells
  i_x1g_a::Int = 1
  i_x1g_b::Int = i_x1_a - 1
  i_x1g_c::Int = i_x1_b + 1
  i_x1g_d::Int = i_x1_b + n_g
  # x2 inside domain
  i_x2_a::Int = 1 + n_g
  i_x2_b::Int = n_x2 + n_g
  # x2 ghost cells
  i_x2g_a::Int = 1
  i_x2g_b::Int = i_x2_a - 1
  i_x2g_c::Int = i_x2_b + 1
  i_x2g_d::Int = i_x2_b + n_g
  # Quantities on the Grid (last axis indices)
  i_rho::Int = 1
  # Conserved variables
  i_rhovx1::Int = 2
  i_rhovx2::Int = 3
  i_rhovx3::Int = 4
  i_rhoE::Int = 5
  # Primitive variables
  i_vx1::Int = 2
  i_vx2::Int = 3
  i_vx3::Int = 4
  i_p::Int = 5
  # Indices normal to interfaces
  i_un::Int = 2
  i_ut::Int = 3
  i_rhoun::Int = 2
  i_rhout::Int = 3

  # Views for convenience
  # Cell-centres
  rho::SubArray{Float64} = @view cc_prim[:, :, i_rho]
  v_x1::SubArray{Float64} = @view cc_prim[:, :, i_vx1]
  v_x2::SubArray{Float64} = @view cc_prim[:, :, i_vx2]
  v_x3::SubArray{Float64} = @view cc_prim[:, :, i_vx3]
  p::SubArray{Float64} = @view cc_prim[:, :, i_p]
  # Left interface
  rho_L::SubArray{Float64} = @view prim_L[:, :, i_rho]
  v_x1_L::SubArray{Float64} = @view prim_L[:, :, i_vx1]
  v_x2_L::SubArray{Float64} = @view prim_L[:, :, i_vx2]
  v_x3_L::SubArray{Float64} = @view prim_L[:, :, i_vx3]
  p_L::SubArray{Float64} = @view prim_L[:, :, i_p]
  # Right interface
  rho_R::SubArray{Float64} = @view prim_R[:, :, i_rho]
  v_x1_R::SubArray{Float64} = @view prim_R[:, :, i_vx1]
  v_x2_R::SubArray{Float64} = @view prim_R[:, :, i_vx2]
  v_x3_R::SubArray{Float64} = @view prim_R[:, :, i_vx3]
  p_R::SubArray{Float64} = @view prim_R[:, :, i_p]

end

"""
TODO:
- Add CRN init for passive scalars
"""