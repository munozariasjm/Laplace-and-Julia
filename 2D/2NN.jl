using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
import ModelingToolkit: Interval, infimum, supremum, Plots
using Pkg

@parameters x y
@variables u(..)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

# 2D PDE
eq  = Dxx(u(x,y)) + Dyy(u(x,y)) ~ 0

# Boundary conditions
bcs = [u(0,y) ~ 10., u(2,y) ~ 0.,
       u(x,0) ~ 0.0, u(x,2) ~ 0.]
# Space and time domains
domains = [x ∈ Interval(0.0,2.0),
           y ∈ Interval(0.0,2.0)]

# Neural network
dim = 2 # number of dimensions
chain = FastChain(FastDense(dim,16,Flux.σ),FastDense(16,16,Flux.σ),FastDense(16,1))
# Initial parameters of Neural network
initθ = Float64.(DiffEqFlux.initial_params(chain))

# Discretization
dx = 0.05
discretization = PhysicsInformedNN(chain,GridTraining(dx),init_params =initθ)

@named pde_system = PDESystem(eq,bcs,domains,[x,y],[u(x, y)])
prob = discretize(pde_system,discretization)

#Optimizer
opt = Optim.BFGS()
losses=[]
#Callback function
cb = function (p,l;ilosses=losses)
    println("Current loss is: $l")
    append!(ilosses,l)
    return false
    
end

res = GalacticOptim.solve(prob, opt, cb = cb, maxiters=1000)
phi = discretization.phi

using Plots

xs,ys = [infimum(d.domain):dx/10:supremum(d.domain) for d in domains]

u_predict = reshape([first(phi([x,y],res.minimizer)) for x in xs for y in ys],(length(xs),length(ys)))

plot(ys, xs, (u_predict), linetype=:contourf,title = "Solución NN")

Plots.plot(losses,label="Pérdida",ylabel="𝕃",xlabel="Iterción entrenamiento")
losses



plot(ys, xs, (u_predict .- Space), linetype=:contourf,title = "Diferencia NN")

using Statistics
mean((u_predict .- Space).^2)
mean((u_predict .- Space)./u_predict)