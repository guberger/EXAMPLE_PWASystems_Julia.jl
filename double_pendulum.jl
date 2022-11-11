using Random
using PyPlot
using DifferentialEquations

Random.seed!(1)

function acrobot!(dx, x, params, ::Any)
    (
        s10, s11, s20, s21, s30, s32, s40, s42,
        a31, a32, a41, a42, b3u, b4u, b3v, b4v
    ) = params
    λ1 = max(s10 + s11*x[1], 0.0)
    λ2 = max(s20 + s21*x[1], 0.0)
    λ3 = max(s30 + s32*x[2], 0.0)
    λ4 = max(s40 + s42*x[2], 0.0)
    dx[1] = x[3]
    dx[2] = x[4]
    dx[3] = a31*x[1] + a32*x[2] + b3u*(λ2 - λ1) + b3v*(λ4 - λ3)
    dx[4] = a41*x[1] + a42*x[2] + b4u*(λ2 - λ1) + b4v*(λ4 - λ3)
end

g = 9.81
m1, l1, d1, k1 = 0.5, 0.5, 0.2, 10
m2, l2, d2, k2 = 0.1, 1, 0.05, 10

a31 = +g/l1
a32 = -(g*m2)/(l1*m1)
a41 = -g/l1
a42 = (g*(l1*m1 + l1*m2 + l2*m2))/(l1*l2*m1)

a = m1*(l1^2) + m2*(l1^2)
b = m2*(l2^2)
c = m2*l1*l2
Ma = [(a + b + 2*c) (b + c); (b + c) b]
iMa = inv(Ma)

b3u = iMa[1, 1]
b4u = iMa[2, 1]
b3v = iMa[1, 2]
b4v = iMa[2, 2]

s10 = -k1*d1
s11 = +k1
s20 = -k1*d1
s21 = -k1
s30 = -k2*d2
s32 = +k2
s40 = -k2*d2
s42 = -k2

params = (
    s10, s11, s20, s21, s30, s32, s40, s42,
    a31, a32, a41, a42, b3u, b4u, b3v, b4v
)

tspan = (0.0, 50.0)
tdom = range(tspan..., length=25)
sols = Vector{Vector{Float64}}[]

for i in 1:10
    x_init = [
        rand()*0.1 - 0.05,
        rand()*0.1 - 0.05,
        rand()*0.4 - 0.2,
        rand()*0.2 - 0.1
    ]
    prob = ODEProblem(acrobot!, x_init, tspan, params)
    sol = solve(prob, saveat=tdom)
    λ1 = (s10 .+ s11*sol[1, :]) .> eps()
    λ2 = (s20 .+ s21*sol[1, :]) .> eps()
    push!(sols, sol[:])
end

fig = figure(figsize=(14, 7))
ax_ = fig.subplots(2, 2)
for i = 1:4
    ax_[i].set_title(string("x[", i, "]"))
    for sol in sols
        ax_[i].plot(tdom, getindex.(sol, i), marker=".")
    end
end

ax_[1].plot(tspan, (-s10/s11, -s10/s11), c="k")
ax_[1].plot(tspan, (-s20/s21, -s20/s21), c="k")
ax_[2].plot(tspan, (-s30/s32, -s30/s32), c="k")
ax_[2].plot(tspan, (-s40/s42, -s40/s42), c="k")

fig = figure()
ax = fig.add_subplot(projection="3d")

for sol in sols
    for x in sol
        λ1 = max(s10 + s11*x[1], 0.0)
        λ2 = max(s20 + s21*x[1], 0.0)
        λ3 = max(s30 + s32*x[2], 0.0)
        λ4 = max(s40 + s42*x[2], 0.0)
        dx3 = a31*x[1] + a32*x[2] + b3u*(λ2 - λ1) + b3v*(λ4 - λ3)
        dx4 = a41*x[1] + a42*x[2] + b4u*(λ2 - λ1) + b4v*(λ4 - λ3)
        ax.plot(x[1], x[2], dx3, marker=".", c="k")
    end
end

IO = open("data/double_pendulum/data.txt", "w")
for sol in sols
    for x in sol
        λ1 = max(s10 + s11*x[1], 0.0)
        λ2 = max(s20 + s21*x[1], 0.0)
        λ3 = max(s30 + s32*x[2], 0.0)
        λ4 = max(s40 + s42*x[2], 0.0)
        dx3 = a31*x[1] + a32*x[2] + b3u*(λ2 - λ1) + b3v*(λ4 - λ3)
        dx4 = a41*x[1] + a42*x[2] + b4u*(λ2 - λ1) + b4v*(λ4 - λ3)
        println(IO, ((x[1], x[2]), (dx3, dx4)))
    end
end
close(IO)