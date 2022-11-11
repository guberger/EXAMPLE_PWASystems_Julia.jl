using Random
using PyPlot
using DifferentialEquations

Random.seed!(1)

function pendulum!(dx, x, params, ::Any)
    s10, s11, s20, s21, a21, b2u = params
    λ1 = max(s10 + s11*x[1], 0.0)
    λ2 = max(s20 + s21*x[1], 0.0)
    dx[1] = x[2]
    dx[2] = a21*x[1] + b2u*(λ2 - λ1) - 0.05*x[2]
end

g, m1, l1, d, k = (9.81, 0.5, 0.5, 0.2, 5)

a21 = +g/l1

Ma = [m1*(l1^2);;]
iMa = inv(Ma)

b2u = iMa[1, 1]

s10 = -k*d
s11 = +k
s20 = -k*d
s21 = -k

params = s10, s11, s20, s21, a21, b2u

tspan = (0.0, 10.0)
tdom = range(tspan..., length=100)
sols = Vector{Vector{Float64}}[]

for i in 1:10
    x_init = [rand()*0.1 - 0.05, 0.0]
    prob = ODEProblem(pendulum!, x_init, tspan, params)
    sol = solve(prob, saveat=tdom)
    push!(sols, sol[:])
end

fig = figure(figsize=(14, 7))
ax_ = fig.subplots(2)
for i = 1:2
    ax_[i].set_title(string("x[", i, "]"))
    for sol in sols
        ax_[i].plot(tdom, getindex.(sol, i), marker=".")
    end
end

ax_[1].plot(tspan, (-s10/s11, -s10/s11), c="k")
ax_[1].plot(tspan, (-s20/s21, -s20/s21), c="k")

fig = figure()
ax = fig.add_subplot()

for sol in sols
    ax.plot(getindex.(sol, 1), getindex.(sol, 2))
end

fig = figure()
ax = fig.add_subplot()

for sol in sols
    for x in sol
        λ1 = max(s10 + s11*x[1], 0.0)
        λ2 = max(s20 + s21*x[1], 0.0)
        dx2 = a21*sin(x[1]) + b2u*(λ2 - λ1)
        ax.plot(x[1], dx2, marker=".", c="k")
    end
end