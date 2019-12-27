# module MyModule

using DifferentialEquations
# using Plots
include("integrate.jl")
using Makie
using AbstractPlotting
using GeometryTypes

print(">>> ")

# plotlyjs()

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

tmax = 5.0
dt = 0.05

sol = Integrate.get_solution(0.0, tmax, dt)

angle = Array{Float64, 1}(undef, size(sol)[1])
position = Array{Float64, 1}(undef, size(sol)[1])
angle_it = 1
pos_it = 1
for i in eachindex(sol)
    global angle_it
    global pos_it
    if (div(i, size(angle)[1]) == 0)
        angle[angle_it] = sol[i]
        angle_it = angle_it + 1
    end
    if (div(i, size(angle)[1]) == 2)
        position[pos_it] = sol[i]
        pos_it = pos_it + 1
    end
end

function convert_angle_to_node(t, v)
    ind = Int(div(t, dt)) + 1
    if (ind >= size(angle)[1])
        return (sin(angle[size(angle)[1]]), cos(angle[size(angle)[1]]))
    else
        return (sin(angle[ind]) + position[ind], cos(angle[ind]))
    end
end

function convert_pos_to_node(t, v)
    ind = Int(div(t, dt)) + 1
    if (ind >= size(position)[1])
        return (position[size(position)[1]], 0.0)
    else
        return (position[ind], 0.0)
    end
end

scene = Scene(resolution = (500, 500))

x = collect(-2:1:2)
y = Array{Float64, 1}(undef, size(x)[1])
for i = 1 : size(y)[1]
    y[i] = 0.0
end

time_node = Node(0.0)
p1 = scatter!(scene, lift(t->convert_angle_to_node.(t, range(0, stop = 1, length = 2)), time_node))[end]
p2 = scatter!(scene, lift(t->convert_pos_to_node.(t, range(0, stop = 1, length = 2)), time_node))[end]

points = lift(p1[1], p2[1]) do pos1, pos2
    map((a, b)-> (a, b), pos1, pos2)
end
linesegments!(scene, points)

lines!(y, x, color = :white)
lines!(x, y, color = :white)

while (true)
    record(scene, "output.mp4", range(0, stop = tmax - dt, length = Int(div(tmax, dt)))) do i
        push!(time_node, i)
    end
end
