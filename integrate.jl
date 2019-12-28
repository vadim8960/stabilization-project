module Integrate

using PyCall

@pyimport scipy.integrate as integrate

struct VariablesODE
    u0::Array{Float64, 1}
    L::Float64
    m::Float64
    M::Float64
end

const g = 9.81
vars = VariablesODE([pi/4.0, 0.0, 0.0, 0.0], 1.0, 1.0, 1.0)

function set_variables(var::VariablesODE)
    global vars
    vars = var
end

function func_stabilization(u, t)
    du = zeros(size(u)[1])
    th = u[1]
    b = vars.m / (vars.m + vars.M)
    up = 50 * th + 15 * u[2] + 3.1 * u[3] + 4.8 * u[4]
    du[1] = u[2]
    du[2] = (g * sin(th) - up * cos(th)) / vars.L
    du[3] = u[4]
    du[4] = up

    return du
end

function func_free_vibrations(u, t)
    du = zeros(size(u)[1])
    th = u[1]
    b = vars.m / (vars.m + vars.M)
    du[1] = u[2]
    du[2] = (g * sin(th) + b * vars.L * du[1]^2 * sin(th) * cos(th)) /
            (vars.L * (1 + b * cos(th)^2)) - du[1]
    du[3] = u[4]
    du[4] = b * (vars.L * du[1]^2 * sin(th) - g * sin(th) * cos(th)) /
            (1 + b * cos(th)^2) - du[3]

    return du
end

function get_solution_stab(t0, tmax, dt)
    return integrate.odeint(func_stabilization, vars.u0, range(t0, stop = tmax - dt, step = dt))
end

function get_solution_free(t0, tmax, dt)
    return integrate.odeint(func_free_vibrations, vars.u0, range(t0, stop = tmax - dt, step = dt))
end

end
