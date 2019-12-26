# module MyModule

using DifferentialEquations
using Plots
include("integrate.jl")
# using Makie

plotlyjs()

# function func(du, u, p, t)
#     th = u[1]
#     du[1] = u[2]
#     du[2] = (g * sin(th) + b * L * du[2]^2 * sin(th) * cos(th)) / (L * (1 + b * cos(th)^2))
#     du[3] = u[4]
#     du[4] = b * (L * du[2]^2 * sin(th) - g * sin(th) * cos(th)) / (1 + b * cos(th)^2)
# end
#
# pror = ODEProblem(func, u0, tspan)
# sol = solve(pror, Tsit5(), reltol = 1e-3)

sol = Integrate.get_solution()

angle = Array{Float64, 1}(undef, div(size(sol)[1], 4))
position = Array{Float64, 1}(undef, div(size(sol)[1], 4))
angle_it = 1
pos_it = 1
for i = 1 : size(sol)[1]
    global angle_it
    global pos_it
    if ((i - 1) % 4 == 0)
        angle[angle_it] = sol[i]
        angle_it = angle_it + 1
    end
    if ((i - 1) % 4 == 2)
        position[pos_it] = sol[i]
        pos_it = pos_it + 1
    end
end

plot(angle)
plot!(position)
