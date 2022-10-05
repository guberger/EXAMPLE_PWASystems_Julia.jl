using Random
using PyPlot
using DifferentialEquations

Random.seed!(0)

function cartpole!(dx, x, params, ::Any)
    g, mp, mc, l, d, k1, k2, K1, K2, K3, K4, L1, L2 = params
    y = x[1] - l*x[2]
    λ1 = max(k1*(+y - d), 0.0)
    λ2 = max(k2*(-y - d), 0.0)
    u = K1*x[1] + K2*x[2] + K3*x[3] + K4*x[4] + L1*λ1 + L2*λ2
    dx[1] = x[3]
    dx[2] = x[4]
    dx[3] = (g*mp/mc)*x[2] + (1/mc)*u
    dx[4] = (g*(mp + mc)/(l*mc))*x[2] + (1/(l*mc))*u + (1/(l*mp))*(λ1 - λ2)
end

g, mp, mc, l, d, k1, k2 = (
    9.81, 0.1, 1, 0.5, 0.1, 10, 10
)
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
params = (g, mp, mc, l, d, k1, k2, K1, K2, K3, K4, L1, L2)

tspan = (0.0, 5.0)
dt = 0.05
tdom = tspan[1]:dt:tspan[2]
sols = Vector{Vector{Float64}}[]

for i in 1:10
    x_init = [
        rand()/5 - 0.1,
        0.0,
        rand()*2 - 1.0,
        rand()*2 - 1.0
    ]
    prob = ODEProblem(cartpole!, x_init, tspan, params)
    sol = solve(prob, saveat=tdom)
    push!(sols, sol[:])
end

fig = figure(figsize=(10, 7))
ax_ = fig.subplots(2, 2)
for i = 1:4
    ax_[i].set_title(string("x[", i, "]"))
    for sol in sols
        ax_[i].plot(tdom, getindex.(sol, i), marker=".")
    end
end

IO = open("data/cartpole/sols.txt", "w")

for sol in sols
    println(IO, sol)
end

close(IO)