#!/usr/bin/env julia

module HydroChem

export hydro, Grid, lax, hll, hllc, euler, RK2, RK3,
  init_sod_x1!, init_sod_x2!, init_ball!, init_KH!, fill_trans_bc, fill_periodic_bc,
  plot_curve, plot_heat, plot_curve_or_heat, plot_heat_four_panels,
  plot_standard_sod
# plot_standard_sod_y,
# hll_vxonly, hllc_vxonly, hllc2, hllc3, hllc3_vy


include("grid.jl")
include("io.jl")
include("flux.jl")
include("bc.jl")
include("init.jl")
include("reconstruction.jl")
include("solver.jl")
# include("hllc3.jl")
include("integrator.jl")
include("plot.jl")
include("jcrn.jl")

const mass_hydrogen = 1.67262171e-24  # [g]
const mm00_abundances = Dict(
  "H" => 12,
  "H2" => -4,
  "OH" => -12,
  "C" => 8.43,  # solar
  "O" => 8.69,  # solar
  "N" => 7.83,
  "CH" => -12,
  "CO" => -12,
  "CN" => -12,
  "NH" => -12,
  "NO" => -12,
  "C2" => -12,
  "O2" => -12,
  "N2" => -12,
  "M" => 11,
  "X" => 1,
  "Y" => 2,
  "Z" => 3
)

""" A general-purpose hydro solver

hydro(dim, nx, tend, folder, init;
      solver::Function=hll,
      order::Int=2,
      integrator::Function=RK3,
      fillbc::Function=fill_trans_bc,
      plotit::Function=plot_curve_or_heat,
      dtout::Float64=0.01,
      ny::Int64=-1,
      storealldata::Bool=false,
      restart=-1,
      islog::Bool=true,
      verbose::Bool=false)

Run a hydrodynamics simulation and save the figures in folder.

Arguments
---------
    dim: (Int) number of dimentions, either 1 or 2
    nx: (Int) number of pixels in the first dimention (use ny= to set the number of pixels in the second dimension)
    tend: (Float) the end time
    folder: (String) the folder to store the outputs (figures, logs, and data)
    init: (Function) a function to setup the initial condition. Choose one of init_sod, init_sod_y, init_ball, init_KH.
    solver: (Function, default: hll) the Riemann solver. Choose one of hll, hllc, lax.
    order: (Int, default: 2) 1 for first order accuracy in space and 2 for second order accuracy.
    integrator: (Function, default: RK3) one of euler, RK3, RK3
    fillbc: (Function) a function to set the boundary condition. One of fill_trans_bc, fill_periodic_bc.
    plotit: (Function) plotting function. One of the following functions.
        plot_curve: plot curves of the density, velocity, pressure, and specific energy
        plot_heat: plot heat map of the density
        plot_heat_four_panels: plot heat maps of the density, x-velocity, y-velocity, and pressure
        plot_curve_or_heat (default): automatically pick plot_curve if dim=1 or plot_heat if dim=2
        plot_standard_sod: on top of plot_curve, plot the analytic solution to the standard Sod shock tube problem (rho = 1.0, 0.125, velocity = 0.0, 0.0, pressure = 1.0, 0.1 for the left and right state)
    dtout: (Float, default: 0.01) the time interval between two outputs.
    ny: (Int, default: -1) the number of pixels in the second dimention. A negative value (default) means ny = nx.
    storealldata: (Bool, default: false) toggle store reloadable data for each outputs set by dtout. If false, will only store the data for the last output, or the last output before you interrupt with Ctrl-c.
    restart: (Int, default: -1) the index of the output from which to restart the simulation. Default value -1 means starting from time 0. 
    islog: (Bool, default: true) toggle creating a log file
    verbose: (Bool, default: false) toggle verbose printing

Examples
--------

    julia> hydro(1, 512, 0.1, "tmp", init_sod)

    julia> hydro(1, 512, 0.2, "tmp", init_sod, restart=10)

    julia> hydro(2, 128, 1.0, "/tmp/Hydro/KH-128", init_KH; dtout=0.005, verbose=true)

"""
function hydro(; n_x1::Int=128, n_x2::Int=1, n_g::Int=2,
  tend::Float64=0.1, folder::String=".", init::Function=init_sod_x1!,
  solver::Function=hll,
  order::Int=2,
  integrator::Function=RK3,
  fillbc::Function=fill_trans_bc,
  plotit::Function=plot_curve_or_heat,
  dtout::Float64=0.01,
  network_file::Union{String,Nothing}=nothing,
  storealldata::Bool=false,
  restart=-1,
  islog::Bool=false,
  verbose::Bool=false)

  msg = """
Running Hydro.jl, a modular 1- and 2-dimensional hydrodynamic
code written in pure Julia.
Author: Chong-Chong He (che1234@umd.edu)

Solving the following problem with parameters:
Initial condition: $(init)
nx = $(n_x1)
dtout = $(dtout)
tend = $(tend)
Riemann solver: $(solver)
Integrator: $(integrator)
Boundary condition: $(fillbc)
"""
  # Plotting function: $(plotit)
  # """
  print(msg)
  @debug "Debug enabled"
  run(`mkdir -p $folder`)
  if islog
    logfile = "$(folder)/log.o"
    if isfile(logfile)
      for i = 1:1000
        logfile = "$(folder)/log.o$(i)"
        if !isfile(logfile)
          break
        end
      end
    end
  end
  if islog
    fw = open(logfile, "w")
    write(fw, msg)
  end
  g = Grid(n_x1=n_x1, n_x2=n_x2, n_g=n_g)
  @info "n_x1 = $(n_x1), n_x2 = $(n_x2), n_g = $(n_g)"
  if order == 1
    reconstruct = reconstruct1st
  elseif order == 2
    reconstruct = reconstruct2nd
  end
  function rebuild(g)
    fillbc(g)               # fill bc for conservatives
    cons2prim!(g)            # convert conservatives to primitives
    return
  end

  datadir = "$(folder)/data"
  if restart < 0
    # New simulation
    g = init(g)         # fill I.C. for prims (without boundary)
    println("Grid initialised")
    # Add chemical species
    if !isnothing(network_file)
      @show size(g.cc_prim)
      g = add_chemical_species!(g, network_file, mm00_abundances)
      @show size(g.cc_prim)
      exit()
    end
    # @show g.cc_prim[g.mid_x1, g.mid_x2, :]
    g = prim2cons!(g)
    println("prim2cons! done")
    # @show g.cc_prim[g.mid_x1, g.mid_x2, :]
    fillbc(g)
    println("fillbc done")
    # @show g.cc_prim[g.mid_x1, g.mid_x2, :]
    g = cons2prim!(g)
    println("cons2prim! done")
    # @show g.cc_prim[g.mid_x1, g.mid_x2, :]
    println("Problem fully initialised.")
    fcount = 0
    run(`mkdir -p $(datadir)`)
    if storealldata
      fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
      grid_write(g, fn)
    end
  else
    # Load from existing file
    # @assert datadir != "none"
    fcount = restart
    fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
    grid_read(g, fn)
  end

  # plot the first snapshot
  println("Initial plot")
  error = plotit(g, fn="$(folder)/hydro_$(lpad(fcount, 5, '0')).png")
  println("Plotted first snapshot.")
  count = 0
  msg = "\ncount = $(count), fcount = $(fcount), t = $(g.t)"
  if plotit == plot_standard_sod
    msg = msg * ", relative error = $(error)"
  end
  print(msg * "\n")
  if islog
    write(fw, msg * "\n")
  end

  # evolve and plot snapshots
  tout = g.t + dtout
  dt = 1.0
  println("Starting evolution from $(g.t) to $(tend) with initial dt = $dt")
  try
    while true
      # v = @view g.u[:, :, 4]
      # println("count = $(count)")
      # println(maximum(v), " ", minimum(v))
      # println()
      dt = calc_dt(g; C=0.8)
      if isnan(dt)
        println("Failed at count = $(count), t = $(g.t): dt = $(dt)")
        # println()
        # v = @view g.u[:, :, 4]
        # println(maximum(v), " ", minimum(v))
        # println()
        # v = g.cs
        # println(maximum(v), " ", minimum(v))
        # println()
        # v = g.pressure
        # println(maximum(v), " ", minimum(v))
        # println()
        # v = g.epsilon
        # println(maximum(v), " ", minimum(v))
        # println()
        # v = g.E
        # println(maximum(v), " ", minimum(v))
        return
      end
      @info "count = $(count), t = $(g.t), dt = $(dt)"
      if verbose
        println("count = $(count), t = $(g.t), dt = $(dt)")
      end
      integrator(g, dt, solver, reconstruct, rebuild)
      count += 1
      if g.t > tout || g.t ≈ tout
        # Final plot, then save data and exit
        fcount += 1
        # plotit(g, fn="$(folder)/f-nx$(nx)-$(lpad(fcount, 5, '0')).png")
        error = plotit(g, fn="$(folder)/hydro_$(lpad(fcount, 5, '0')).png")
        msg = "count = $(count), fcount = $(fcount), " *
              "t = $(g.t), tout = $(tout)"
        if !(error == nothing)
          msg = msg * ", relative error = $(error)"
        end
        print(msg * "\n")
        if islog
          write(fw, msg * "\n")
        end
        if storealldata || (g.t > tend || g.t ≈ tend)
          # write jld data
          fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
          grid_write(g, fn)
        end
        if g.t > tend || g.t ≈ tend
          break
        end
        tout += dtout
      end
    end
  catch err
    if isa(err, InterruptException)
      fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
      msg = "Interrupted, saving $(fn)\n"
      print(msg)
      grid_write(g, fn)
      if islog
        write(fw, msg)
      end
      # return
    else
      rethrow()
    end
  end
  println("Simulation done.")
  if islog
    write(fw, "Simulation done.\n")
    close(fw)
  end
  return 0
end

end

HydroChem.hydro(n_x1=128, folder="test", network_file="../res/solar_co_w05.ntw")

"""
Current issues:
  - Change arrays to be arr[quantity, y_idx, x_idx] to exploit column-major
    operations
  - Views in Grid don't seem to update the original array in Grid!
"""