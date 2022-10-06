""" Set the 2d initial conditions: a 2D version of the 1D sod shocktube """
function init_sod_x1!(g::Grid)
    println("Inside init_sod_x1!")
    rho1 = 1.0
    rho2 = 0.125
    p1 = 1.0
    p2 = 0.1

    # @. g.rho[1:g.mid_x1, :] = rho1
    # @. g.p[1:g.mid_x1, :] = p1
    # @. g.rho[g.mid_x1+1:end, :] = rho2
    # @. g.p[g.mid_x1+1:end, :] = p2
    # @. g.v_x1 = 0.0
    # @. g.v_x2 = 0.0
    @. g.cc_prim[1:g.mid_x1, :, g.i_rho] = rho1
    @. g.cc_prim[1:g.mid_x1, :, g.i_p] = p1
    @. g.cc_prim[g.mid_x1+1:end, :, g.i_rho] = rho2
    @. g.cc_prim[g.mid_x1+1:end, :, g.i_p] = p2
    # @. g.v_x1 = 0.0
    # @. g.v_x2 = 0.0
    return g
end


""" Set the 2d initial conditions: a 2D version of the 1D sod shocktube in the y dimension """
function init_sod_x2!(g::Grid)
    mid = floor(Int64, g.x2_len / 2)
    rho1 = 1.0
    rho2 = 0.125
    p1 = 1.0
    p2 = 0.1

    g.rho[:, 1:mid] .= rho1
    g.p[:, 1:mid] .= p1
    g.rho[:, (mid+1):end] .= rho2
    g.p[:, (mid+1):end] .= p2
    g.v_x1 .= 0.0
    g.v_x2 .= 0.0
    return g
end


""" Set the 2d sod shock tube: a ball of radius 0.2 at center with
density 1.0, pressure 1.0. Outside has density 0.125 and pressure 0.1 """
function init_ball!(g::Grid)
    radius = 0.2
    midx = floor(Int64, g.x1_len / 2)
    midy = floor(Int64, g.x2_len / 2)
    rho1 = 1.0
    rho2 = 0.125
    p1 = 1.0
    p2 = 0.1

    cx = g.x1[midx]
    cy = g.x2[midy]
    for i = 1:g.x1_len, j = 1:g.x2_len
        if (g.x1[i] - cx)^2 + (g.x2[j] - cy)^2 <= radius^2
            g.rho[i, j] = rho1
            g.p[i, j] = p1
        else
            g.rho[i, j] = rho2
            g.p[i, j] = p2
        end
    end
    g.v_x1 .= 0.0
    g.v_x2 .= 0.0
    return g
end


""" Kelvin-Helmholtz instability """
function init_KH!(g::Grid)
    rho1 = 1.0
    rho2 = 2.0
    v1 = -0.5
    v2 = 0.5
    p = 2.5
    Amp = 0.01

    for j = 1:g.x2_len
        if abs(g.x2[j] - 0.5) > 0.25
            g.rho[:, j] .= rho1
            g.v_x1[:, j] .= v1
        else
            g.rho[:, j] .= rho2
            g.v_x1[:, j] .= v2
        end
    end
    g.p .= p
    # # initial perturbation, old version
    # g.u[:, :, 2] .+= 0.01 * (1 .- 2 * rand(g.x1_len, g.x2_len)) .* g.u[:, :, 1]
    # g.u[:, :, 3] .+= 0.01 * (1 .- 2 * rand(g.x1_len, g.x2_len)) .* g.u[:, :, 1]
    # initial perturbation: sin wave, Springel 2009
    σ = 0.05
    for j = g.i_x2_a:g.i_x2_b
        y = g.x2[j]
        # vy = 0.01 * sin.(4 * pi * collect(1:g.x1_len) / g.x1_len)
        # expo = - (y - 0.25)^2 / (2 * σ^2) + (y - 0.75)^2 / (2 * σ^2)
        g.v_x2[:, j] .= 0.1 * sin.(4 * pi * g.x) .* (exp(-(y - 0.25)^2 / (2 * σ^2)) .+ exp.(-(y - 0.75)^2 / (2 * σ^2)))
    end
    return g
end


""" Kelvin-Helmholtz instability """
function init_KH2(g::Grid)
    rho1 = 1.0
    rho2 = 2.0
    v1 = 0.5
    v2 = -0.5
    p = 2.5
    Amp = 0.01

    for j = 1:g.x2_len
        if abs(g.x2[j] - 0.5) > 0.25
            g.rho[:, j] .= rho1
            g.v_x1[:, j] .= v1
        else
            g.rho[:, j] .= rho2
            g.v_x1[:, j] .= v2
        end
    end
    Lx = g.x1[g.i_x1_b] - g.x1[g.i_x1_a]
    Ly = g.x2[g.i_x2_b] - g.x2[g.i_x2_a]
    for j = g.i_x2_a:g.i_x2_b, i = g.i_x1_a:g.i_x1_b
        g.v_x1[i, j] += Amp * (1.0 + cos(2pi * 2 * g.x1[i] / Lx) * (1.0 + cos(2pi * 2 * g.x2[j] / Ly)))
    end
    g.v_x2 .= 0.0
    g.p .= p

    return
end


# random velocities with amplitude 0.01
function init_KH_rand(g::Grid)
    rho1 = 1.0
    rho2 = 2.0
    v1 = -0.5
    v2 = 0.5
    p = 2.5
    Amp = 0.01

    for j = 1:g.x2_len
        if abs(g.x2[j] - 0.5) > 0.25
            g.rho[:, j] .= rho1
            g.v_x1[:, j] .= 0.1
        else
            g.rho[:, j] .= rho2
            g.v_x1[:, j] .= 0.1
        end
    end
    g.v_x1 .+= Amp .* (2 * rand(size(g.v_x1)) - 1.0)
    g.v_x2 .= Amp .* (2 * rand(size(g.v_x1)) - 1.0)
    g.p .= p
end
