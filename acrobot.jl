using Random
using PyPlot
using DifferentialEquations

Random.seed!(1)

function acrobot!(dx, x, params, ::Any)
    (
        s10, s11, s20, s21,
        a31, a32, a41, a42, b3, b4, d31, d32, d41, d42,
        K1, K2, K3, K4, L1, L2
    ) = params
    λ1 = max(s10 + s11*x[1], 0.0)
    λ2 = max(s20 + s21*x[1], 0.0)
    u = K1*x[1] + K2*x[2] + K3*x[3] + K4*x[4] + L1*λ1 + L2*λ2
    dx[1] = x[3]
    dx[2] = x[4]
    dx[3] = a31*x[1] + a32*x[2] + b3*u + d31*λ1 + d32*λ2
    dx[4] = a41*x[1] + a42*x[2] + b4*u + d41*λ1 + d42*λ2
end

g, m1, m2, l1, l2, d, k = (9.81, 0.5, 1, 0.5, 1, 1, 1)

a31 = +g/l1
a32 = -(g*m2)/(l1*m1)
a41 = -g/l1
a42 = (g*(l1*m1 + l1*m2 + l2*m2))/(l1*l2*m1)

b3 = -(l1 + l2)/(l2*m1*l1^2)
b4 = ((m1*l1^2) + m2*(l1 + l2)^2)/((l1^2)*(l2^2)*m1*m2)

a = m1*(l1^2) + m2*(l1^2)
b = m2*(l2^2)
c = m2*l1*l2
Ma = [(a + b + 2*c) (b + c); (b + c) b]
iMa = inv(Ma)

d31 = -iMa[1, 1]
d32 = +iMa[1, 1]
d41 = -iMa[2, 1]
d42 = +iMa[2, 1]

s10 = -k*d
s11 = +k
s20 = -k*d
s21 = -k

K1, K2, K3, K4 = (
    73.06931508,
    38.1070114,
    30.4055215,
    18.94774763
)
L1, L2 = (
    -4.13491252,
    4.13495174
)

params = (
    s10, s11, s20, s21,
    a31, a32, a41, a42, b3, b4, d31, d32, d41, d42,
    K1, K2, K3, K4, L1, L2
)

tspan = (0.0, 3.0)
dt = 0.03
tdom = tspan[1]:dt:tspan[2]
sols = Vector{Vector{Float64}}[]
contacts = Vector{Vector{Float64}}[]

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
    push!(contacts, map(collect, zip(λ1, λ2)))
end

fig = figure(figsize=(14, 7))
ax_ = fig.subplots(2, 3)
for i = 1:4
    ax_[i].set_title(string("x[", i, "]"))
    for sol in sols
        ax_[i].plot(tdom, getindex.(sol, i), marker=".")
    end
end
for i = 1:2
    ax_[4 + i].set_title(string("λ", i))
    for contact in contacts
        ax_[4 + i].plot(tdom, getindex.(contact, i), marker=".")
    end
end

ax_[1].plot(tspan, (-s10/s11, -s10/s11), c="k")
ax_[1].plot(tspan, (-s20/s21, -s20/s21), c="k")

IO = open("data/acrobot/sols.txt", "w")
for sol in sols
    println(IO, sol)
end
close(IO)

IO = open("data/acrobot/contacts.txt", "w")
for contact in contacts
    println(IO, contact)
end
close(IO)