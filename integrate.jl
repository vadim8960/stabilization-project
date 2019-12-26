module Integrate

using PyCall

py"""

import scipy.integrate as int
import numpy as np
from math import *

u0 = np.array([pi/4.0, 0.0, 0.0, 0.0])
g = 9.8
L = 1
m = 1
M = 1
b = m / (m + M)

def func(u, t):
    du = np.zeros_like(u)
    th = u[0]
    du[0] = u[1]
    du[1] = (g * sin(th) + b * L * du[0]**2 * sin(th) * cos(th)) / (L * (1 + b * cos(th)**2))
    du[2] = u[3]
    du[3] = b * (L * du[0]**2 * sin(th) - g * sin(th) * cos(th)) / (1 + b * cos(th)**2)

    return du

def get_solution():
    anw = int.odeint(func, u0, np.arange(0, 10, 0.05))
    return anw.reshape(len(anw) * len(anw[1]))

"""

function get_solution()
    py"get_solution"()
end

end
