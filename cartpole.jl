using Random
using PyPlot
using DifferentialEquations

Random.seed!(1)

function cartpole!(dx, x, params, ::Any)
    (
        s10, s11, s12, s20, s21, s22,
        a32, a42, b3, b4, d41, d42,
        K1, K2, K3, K4, L1, L2
    ) = params
    λ1 = max(s10 + s11*x[1] + s12*x[2], 0.0)
    λ2 = max(s20 + s21*x[1] + s22*x[2], 0.0)
    u = K1*x[1] + K2*x[2] + K3*x[3] + K4*x[4] + L1*λ1 + L2*λ2
    dx[1] = x[3]
    dx[2] = x[4]
    dx[3] = a32*x[2] + b3*u
    dx[4] = a42*x[2] + b4*u + d41*λ1 + d42*λ2
end

g, mp, mc, l, d, k1, k2 = (9.81, 0.1, 1, 0.5, 0.1, 10, 10)

a32 = g*mp/mc
a42 = g*(mp + mc)/(l*mc)
b3 = 1/mc
b4 = 1/(l*mc)
d41 = +1/(l*mp)
d42 = -1/(l*mp)

s10 = -k1*d
s11 = +k1
s12 = -k1*l
s20 = -k2*d
s21 = -k2
s22 = +k2*l

K1, K2, K3, K4 = (
    1.600433659949915,
    -39.029745978214010,
    1.884785337180030,
    -4.987108356077495
)
L1, L2 = (
    -12.909522105700301,
    12.913673613969927
)

params = (
    s10, s11, s12, s20, s21, s22,
    a32, a42, b3, b4, d41, d42,
    K1, K2, K3, K4, L1, L2
)

tspan = (0.0, 5.0)
dt = 0.05
tdom = tspan[1]:dt:tspan[2]
sols = Vector{Vector{Float64}}[]
contacts = Vector{Vector{Float64}}[]

for i in 1:10
    x_init = [
        rand()/5 - 0.1,
        0.0,
        rand()*2 - 1.0,
        rand()*2 - 1.0
    ]
    prob = ODEProblem(cartpole!, x_init, tspan, params)
    sol = solve(prob, saveat=tdom)
    λ1 = (s10 .+ s11*sol[1, :] + s12*sol[2, :]) .> eps()
    λ2 = (s20 .+ s21*sol[1, :] + s22*sol[2, :]) .> eps()
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

IO = open("data/cartpole/sols.txt", "w")
for sol in sols
    println(IO, sol)
end
close(IO)

IO = open("data/cartpole/contacts.txt", "w")
for contact in contacts
    println(IO, contact)
end
close(IO)