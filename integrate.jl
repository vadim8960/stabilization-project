module Integrate

using PyCall

py"""

import scipy.integrate as int
import numpy as np
from math import *

u0 = np.array([pi/10, 0.0, 0.0, 0.0])
g = 9.8
L = 1
m = 0.1
M = 10
b = m / (m + M)

def func(u, t):
    du = np.zeros_like(u)
    th = u[0]
    # du[0] = u[1]
    # du[1] = (g * sin(th) + b * L * du[0]**2 * sin(th) * cos(th)) / (L * (1 + b * cos(th)**2)) - du[0]
    # du[2] = u[3]
    # du[3] = b * (L * du[0]**2 * sin(th) - g * sin(th) * cos(th)) / (1 + b * cos(th)**2) - du[2]
    up = 49 * th + 15 * u[1] + 3.1 * u[2] + 4.8 * u[3]
    du[0] = u[1]
    du[1] = (g * sin(th) - up * cos(th)) / L
    du[2] = u[3]
    du[3] = up
	# du[0] = u[1]
	# du[1] = (g * sin(th) - 1 * cos(th)) / L
	# du[2] = u[3]
	# du[3] = 1

    return du

def get_solution(t0, tmax, dt):
    anw = int.odeint(func, u0, np.arange(t0, tmax, dt))
    return anw.reshape(len(anw) * len(anw[1]))

"""

function get_solution(t0, tmax, dt)
    py"get_solution"(t0, tmax, dt)
end

end
