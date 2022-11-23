using Parameters

include("jcrn.jl")

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

  # Passive scalars
  ps::SubArray{Float64} = @view cc_prim[:, :, 6:nvars]

  # Chemistry (initialised with functions afterwards)
  # TODO: Make a mutable struct to hold this and access that
  network_file::Union{String,Nothing} = nothing
  reaction_network::Union{ReactionSystem,Nothing} = nothing
  species::Union{Array,Nothing} = nothing
  chem_ode_prob::Union{ODEProblem,Nothing} = nothing
end

function resize_arr(g::Grid, old_arr::Array{Float64,3}, nvars_new::Int)
  # Convenience method to create a new array with last dim size nvars_new > 5
  # Used to add passive scalars (initialised to zero) to all relevant arrays,
  # but keeping values from previous Grid for all other variables
  @info "Resizing array..."
  if nvars_new < 5
    println("Error: nvars_new must be greater than 5.")
    exit()
  end
  new_arr::Array{Float64,3} = zeros(Float64, g.x1_len, g.x2_len, nvars_new)
  @. new_arr[:, :, 1:5] = old_arr[:, :, 1:5]
  return new_arr
end

function remake_grid!(g::Grid, n_x1::Int, n_x2::Int, n_g::Int, n_ps::Int)
  # As of current implementation, only has functionality for adding passive
  # scalars to arrays; future functionality should allow for extending and
  # truncating x1, x2 sizes, but this is not handled
  # For passive scalars, overwrite previous grid passive scalars, copying all
  # conserved and primitive variables and allocating space for new scalars
  # Since all passive scalar values are set to zero, the user needs to
  # populate these separately and then recompute left & right interfaces
  if n_x1 != g.n_x1
    println("x1-cells: $(g.n_x1) -> $(n_x1): Currently not supported.")
  end
  if n_x2 != g.n_x2
    println("x2-cells: $(g.n_x2) -> $(n_x2): Currently not supported.")
  end
  if n_g != g.n_g
    println("Ghost cells: $(g.n_g) -> $(n_g): Currently not supported.")
  end
  if n_ps != g.n_ps
    println("Passive scalars: $(g.n_ps) -> $(n_ps): Overwriting previous passive scalars.")
  end
  nvars_new::Int = 5 + n_ps
  g.n_ps = n_ps
  g.nvars = nvars_new
  # Fill conserved and primitive variables from previous grid
  g.cc_cons = resize_arr(g, g.cc_cons, nvars_new)
  g.cc_prim = resize_arr(g, g.cc_prim, nvars_new)

  # Left and right interface
  g.cons_L = resize_arr(g, g.cons_L, nvars_new)
  g.cons_R = resize_arr(g, g.cons_R, nvars_new)
  g.prim_L = resize_arr(g, g.prim_L, nvars_new)
  g.prim_R = resize_arr(g, g.prim_R, nvars_new)

  # Fluxes
  g.f_cons_x1 = resize_arr(g, g.f_cons_x1, nvars_new)
  g.f_cons_x2 = resize_arr(g, g.f_cons_x2, nvars_new)

  g.lu = resize_arr(g, g.lu, nvars_new)
  g.fhll = resize_arr(g, g.fhll, nvars_new)

  # Views
  # Cell-centres
  g.rho::SubArray{Float64} = @view g.cc_prim[:, :, g.i_rho]
  g.v_x1::SubArray{Float64} = @view g.cc_prim[:, :, g.i_vx1]
  g.v_x2::SubArray{Float64} = @view g.cc_prim[:, :, g.i_vx2]
  g.v_x3::SubArray{Float64} = @view g.cc_prim[:, :, g.i_vx3]
  g.p::SubArray{Float64} = @view g.cc_prim[:, :, g.i_p]
  g.ps::SubArray{Float64} = @view g.cc_cons[:, :, 6:g.nvars]
  # Left interface
  g.rho_L::SubArray{Float64} = @view g.prim_L[:, :, g.i_rho]
  g.v_x1_L::SubArray{Float64} = @view g.prim_L[:, :, g.i_vx1]
  g.v_x2_L::SubArray{Float64} = @view g.prim_L[:, :, g.i_vx2]
  g.v_x3_L::SubArray{Float64} = @view g.prim_L[:, :, g.i_vx3]
  g.p_L::SubArray{Float64} = @view g.prim_L[:, :, g.i_p]
  # Right interface
  g.rho_R::SubArray{Float64} = @view g.prim_R[:, :, g.i_rho]
  g.v_x1_R::SubArray{Float64} = @view g.prim_R[:, :, g.i_vx1]
  g.v_x2_R::SubArray{Float64} = @view g.prim_R[:, :, g.i_vx2]
  g.v_x3_R::SubArray{Float64} = @view g.prim_R[:, :, g.i_vx3]
  g.p_R::SubArray{Float64} = @view g.prim_R[:, :, g.i_p]

  # @show any(isnan, g.rho)
  # @show any(isnan, g.p)
  return g
end

function add_chemical_species!(g::Grid, network_file::String, abundances::Dict)
  # Read chemical network file, create reaction system
  # Add chemical species as passive scalars in alphabetical order
  @info "Reading reaction network..."
  g.reaction_network = read_network_file(network_file)
  g.species = str_replace.(species(g.reaction_network))
  sort!(g.species)
  @info "Adding chemical species to Grid as passive scalars"
  g = remake_grid!(g, g.n_x1, g.n_x2, g.n_g, length(g.species))

  # Assign number densities
  @info "Assigning number densities..."
  for i = 1:g.x1_len, j = 1:g.x2_len
    g.ps[i, j, :] = calculate_number_densities_arr(g.rho[i, j], abundances, g.species)
  end
  # Verbose
  for (k, s) in enumerate(g.species)
    @debug "$(s): $(minimum(log10.(g.ps[:, :, k]))) $(maximum(log10.(g.ps[:, :, k])))"
  end
  return g
end
