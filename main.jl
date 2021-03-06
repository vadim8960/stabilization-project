using DifferentialEquations
include("integrate.jl")
using Makie

#                   Variable note                   u0             L    m    M
start_variables = Integrate.VariablesODE([pi/4.0, 0.0, 0.0, 0.0], 1.0, 1.0, 1.0)

Integrate.set_variables(start_variables)

x_axis_min = -(start_variables.L + 5)
x_axis_max = (start_variables.L + 5)

tmax = 10.0
dt = 0.05

sol = vec(Integrate.get_solution_stab(0.0, tmax, dt))

size_meas = Int(div(size(sol)[1], 4))

angle = sol[filter(x->(Int(div(x - 1, size_meas)) == 0), range(1, stop = size(sol)[1], step = 1))]
position = sol[filter(x->(Int(div(x - 1, size_meas)) == 2), range(1, stop = size(sol)[1], step = 1))]

function convert_angle_to_node(t, v)
    ind = Int(div(t, dt)) + 1
    return (start_variables.L * sin(angle[ind]) + position[ind], start_variables.L * cos(angle[ind]))
end

function convert_pos_to_node(t, v)
    ind = Int(div(t, dt)) + 1
    return (position[ind], 0.0)
end

plot_scene = plot(range(0, stop = tmax - dt, step = dt), angle, color = :red)
plot!(plot_scene, range(0, stop = tmax - dt, step = dt), position, color = :blue)

scene = Scene(resolution = (500, 500))

x = collect(x_axis_min:1:x_axis_max)
y = zeros(size(x)[1])

lines!(scene, y, x, color = :lightgray)
lines!(scene, x, y, color = :lightgray)

time_node = Node(0.0)
p1 = scatter!(scene, lift(t->convert_angle_to_node.(t, range(0, stop = 1, length = 2)), time_node))[end]
p2 = scatter!(scene, lift(t->convert_pos_to_node.(t, range(0, stop = 1, length = 2)), time_node))[end]

points = lift(p1[1], p2[1]) do pos1, pos2
    map((a, b)-> (a, b), pos1, pos2)
end
linesegments!(scene, points)

while (true)
    record(vbox(scene, plot_scene), "output.mp4", range(0, stop = tmax - dt, step = dt)) do i
        push!(time_node, i)
    end
end
