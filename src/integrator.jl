""" Time integrator """

function calc_dt_new_test(g::Grid, dt_init)
    C = 0.9
    dt = Inf
    for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
        dt = min(dt, g.dx1 / max(abs(g.v_x1[i, j] + g.cs[i, j]), abs(g.v_x1[i, j] - g.cs[i, j])))
        dt = min(dt, g.dx2 / max(abs(g.v_x2[i, j] + g.cs[i, j]), abs(g.v_x2[i, j] - g.cs[i, j])))
    end
    return min(C * dt, 1.05 * dt_init)
end


function calc_dt(g::Grid; C=0.8)
    themax = 0.0
    for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
        vx = g.v_x1[i, j]
        vy = g.v_x2[i, j]
        cs = g.cs[i, j]
        themax = max(themax,
            max(abs(vx + cs), abs(vx - cs)) / g.dx1 + max(abs(vy + cs), abs(vy - cs)) / g.dx2)
    end
    return C / themax
    # return C / maximum(max.(abs.(g.v_x1 .+ g.cs), abs.(g.v_x1 .- g.cs)) ./ g.dx1 .+ max.(abs.(g.v_x2 .+ g.cs), abs.(g.v_x2 .- g.cs)) ./ g.dx2)
end

# General Euler, 1st order
function euler(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    lu = solver(g, reconstruct)
    @. g.cc_cons = g.cc_cons + dt * lu
    rebuild(g)
    g.t += dt
    return
end

# RK2, general
function RK2(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    uold = copy(g.cc_cons)
    lu = solver(g, reconstruct)
    k1 = dt .* lu
    @. g.cc_cons = g.cc_cons + 0.5 * k1
    rebuild(g)
    lu = solver(g, reconstruct)
    k2 = dt .* lu
    @. g.cc_cons = uold + k2
    rebuild(g)
    g.t += dt
end

# RK2, general
function RK2new(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    uold = copy(g.cc_cons)
    lu = solver(g, reconstruct)
    k1 = dt .* lu
    @. g.cc_cons = g.cc_cons + k1
    rebuild(g)
    lu = solver(g, reconstruct)
    k2 = dt .* lu
    @. g.cc_cons = uold + (k1 + k2) / 2
    rebuild(g)
    g.t += dt
end

# RK3, general
function RK3(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    println("Inside RK3")
    uold = copy(g.cc_cons)
    mat = [1.0 0.0 1.0; 0.75 0.25 0.25; 1/3 2/3 2/3]
    for i = 1:size(mat, 1)
        lu = solver(g, reconstruct)
        # v = @view lu[:, :, 4]
        # println("lu[:, :, 4] = ", maximum(v), " ", minimum(v))
        @. g.cc_cons = mat[i, 1] * uold + mat[i, 2] * g.cc_cons + mat[i, 3] * dt * lu
        rebuild(g)
    end
    g.t += dt
end

"""
TODO:
- Embed JCRN calls after each hydro time step
"""