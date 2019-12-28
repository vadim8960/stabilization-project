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

x_axis_min = -3
x_axis_max = 3

tmax = 10.0
dt = 0.05

sol = vec(Integrate.get_solution_stab(0.0, tmax, dt))

size_meas = Int(div(size(sol)[1], 4))

angle = sol[filter(x->(Int(div(x - 1, size_meas)) == 0), range(1, stop = size(sol)[1], step = 1))]
position = sol[filter(x->(Int(div(x - 1, size_meas)) == 2), range(1, stop = size(sol)[1], step = 1))]

println(angle)

function convert_angle_to_node(t, v)
    ind = Int(div(t, dt)) + 1
    return (sin(angle[ind]) + position[ind], cos(angle[ind]))
end

function convert_pos_to_node(t, v)
    ind = Int(div(t, dt)) + 1
    return (position[ind], 0.0)
end

println(size(angle))
println(size(position))

scene = Scene(resolution = (500, 500))

plot_scene = plot(range(0, stop = tmax - dt, step = dt), angle, color = :red)
plot!(plot_scene, range(0, stop = tmax - dt, step = dt), position, color = :blue)

x = collect(x_axis_min:1:x_axis_max)
y = zeros(size(x)[1])

time_node = Node(0.0)
p1 = scatter!(scene, lift(t->convert_angle_to_node.(t, range(0, stop = 1, length = 2)), time_node))[end]
p2 = scatter!(scene, lift(t->convert_pos_to_node.(t, range(0, stop = 1, length = 2)), time_node))[end]

points = lift(p1[1], p2[1]) do pos1, pos2
    map((a, b)-> (a, b), pos1, pos2)
end
linesegments!(scene, points)

lines!(scene, y, x, color = :white)
lines!(scene, x, y, color = :white)

while (true)
    record(vbox(scene, plot_scene), "output.mp4", range(0, stop = tmax - dt, length = Int(div(tmax, dt)))) do i
        push!(time_node, i)
    end
end
