
###############################################
# Solve laplace equation

# by: Jose Miguel Muñoz
# please contact me if this doesn't works
###############################################

using Plots, Statistics, BenchmarkTools
###############################################
#       BOUNDARY
###############################################
N = 100
a , b = 0,2

Rx = a:(b-a)/(N-1):b
Ry = Rx

boundary_X_Up = map(rx -> 0,Rx)

boundary_X_Down = map(rx -> 0, Rx)

boundary_Y_R = map(rx -> 10,Rx)

boundary_Y_L = map(rx -> 0,Rx)


V = zeros(N, N)

V[1,:] = boundary_X_Up
V[:,1] = boundary_Y_R
V[N,:] = boundary_X_Down
V[:,N] = boundary_Y_L

###############################################
#       SOLVER
###############################################
function solver!(V,n_inter,N=N, tol=1)
    for rep = 1:n_inter
        
        for row = 2:N-1
            for column = 2:N-1
                V[row, column] = mean([V[row-1, column],V[row+1, column],V[row, column-1],V[row, column+1]])
                #if (row+10)^2+(column+10)^2<1000
                #    V[row, column]=1
                #end
    end end
 end
    return V
end

###############################################
#       Ask for the boundaries
###############################################
@benchmark solver!(V,1000)
V = solver!(V,100000)

plot(Ry, Ry, V, linetype=:contourf,title = "Solucion por Gauss_Seidel")

#gr()
#anim = @animate for n = 1:100:1000
#    heatmap(vss[:,:,n])
#end
#gif(anim)
using Statistics
plot(Ry, Ry, (V-Space*1.01), linetype=:contourf,title = "Solucion por Jacobi")
mean(V-Space)
mean(V.^2-Space.^2)
std(V.^2-Space.^2)


#######################################################
#       Convergence Analysis
#######################################################
listconvGauss=[]
for n_reps=500:100:5000
    V = rand(N, N)

    V[1,:] = boundary_X_Up
    V[:,1] = boundary_Y_R
    V[N,:] = boundary_X_Down
    V[:,N] = boundary_Y_L

    V = solver!(V,n_reps)
    append!(listconvGauss,abs(mean(V-Space)))
end
listconvGauss


alistconvSOR=map(x->abs(x),listconvSOR)
alistconvJacobi=map(x->abs(x),listconvJacobi)
alistconvGauss=map(x->abs(x),listconvGauss)
plot()
plot!(Vector(500:100:5000),alistconvGauss,label="Metodo Gauss_Seidel")
plot!(ps,alistconvSOR[begin+6:end],label="Metodo SOR")
plot!(ps,alistconvJacobi.+randn(),label="Metodo Jacobi")
xlabel!("Numero de iteraciones")
ylabel!("Δφ")
ps
alistconvSOR[begin:end-6]

alistconvJacobi
ps