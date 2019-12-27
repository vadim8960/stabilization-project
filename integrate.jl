module Integrate

using PyCall

@pyimport scipy.integrate as integrate

struct VariablesODE
    u0::Array{Float64, 1}
    L::Float64
    m::Float64
    M::Float64
end

vars = VariablesODE([pi/4.0, 0.0], 1.0, 1.0, 1.0)

function f(u0, t)
    du = zeros(Float64, 2)

    du[1] = u0[2]
    du[2] = u0[1] * vars.L

    return du
end

println(integrate.odeint(f, vars.u0, range(0, stop = 1, step = 0.1)))

py"""

import scipy.integrate as int
import numpy as np
from math import *

g = 9.8

u0 = np.array([pi/4.0, 0.0, 0.0, 0.0])
L = 1
m = 0.1
M = 1

def set_phis_params(new_L, new_m, new_M):
    L = new_L
    m = new_m
    M = new_M

def func_stabilization(u, t):
    du = np.zeros_like(u)
    th = u[0]
    b = m / (m + M)
    up = 50 * th + 15 * u[1] + 3.1 * u[2] + 4.8 * u[3]
    du[0] = u[1]
    du[1] = (g * sin(th) - up * cos(th)) / L
    du[2] = u[3]
    du[3] = up

    return du

def func_free_vibrations(u, t):
    du = np.zeros_like(u)
    th = u[0]
    b = m / (m + M)
    du[0] = u[1]
    du[1] = (g * sin(th) + b * L * du[0]**2 * sin(th) * cos(th)) / (L * (1 + b * cos(th)**2)) - du[0]
    du[2] = u[3]
    du[3] = b * (L * du[0]**2 * sin(th) - g * sin(th) * cos(th)) / (1 + b * cos(th)**2) - du[2]

    return du

def get_solution(t0, tmax, dt):
    anw = int.odeint(func_stabilization, u0, np.arange(t0, tmax, dt))
    return anw.reshape(len(anw) * len(anw[1]))

"""

function get_solution(t0, tmax, dt)
    py"get_solution"(t0, tmax, dt)
end

end
