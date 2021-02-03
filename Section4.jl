using DifferentialEquations
using Plots
gr()
theme(:ggplot2)

function ode_func!(du,u,params,t)
    q, b = u
    α, β, δ, λ, v, p, E, a, K, G = params
    du[1] = dq = λ*(1-q)*b*(1-G(v-p+E-a+K)) - β*q # u[1]
    du[2] = db = α*λ*q*(1-b)*(1-G(v-p)) - δ*b # u[2]
end

a = 2; K=1; E=0.5; α=0.5; β=0.5; δ=0.9; v=3; λ=1; p=1000
G(x) = 0 ≤ x ≤ 10 ? 1/10 : 0 # uniform distribution on [0,10]
params = [α, β, δ, λ, v, p, E, a, K, G]
q0=0.4; b0=0.2; u0 = [q0, b0]
tspan=(0.0,10.0)

ode_prob = ODEProblem(ode_func!,u0,tspan,params)

ode_sol = solve(ode_prob)

xyzt = plot(ode_sol, plotdensity=10000,lw=1.5)
xy = plot(ode_sol, plotdensity=10000, vars=(1), label = "q")
xz = plot(ode_sol, plotdensity=10000, vars=(2), label = "b")
xyz = plot(ode_sol, plotdensity=10000, vars=(1,2))
plot(plot(xyzt,xyz),plot(xy, xz, layout=(1,2),w=1), layout=(2,1))


