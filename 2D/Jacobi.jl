###############################################
# Solve laplace equation

# by: Jose Miguel Muñoz
# please contact me if this doesn't works
###############################################

using Plots, Statistics
###############################################
#       BOUNDARY
###############################################
N = 100
a , b = 0, 2

Rx = a:(b-a)/(N-1):b
Ry = a:(b-a)/(N-1):b

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
        nV =  (V)
        for row = 2:N-1
            for column = 2:N-1
                V[row, column] = mean([nV[row-1, column],nV[row+1, column],nV[row, column-1],nV[row, column+1]])
                #if (row+10)^2+(column+10)^2<1000
                #    V[row, column]=1
                #end
    end end
    V =  (nV)
 end
    return V
end

###############################################
#       Ask for the boundaries
###############################################
V = solver!(V,100000)
@benchmark solver!(V,1000)
#contour(V,levels=100)
plot(Ry, Ry, V, linetype=:contourf,title = "Solucion por Jacobi")

#gr()
#anim = @animate for n = 1:100:1000
#    heatmap(vss[:,:,n])
#end
#gif(anim)
##############################
#   Errores
##############################
using Statistics
plot(Ry, Ry, (V-Space*1.01), linetype=:contourf,title = "Solucion por Jacobi")

std(V.^2-Space.^2)


#######################################################
#       Convergence Analysis
#######################################################
listconvJacobi=[]
for n_reps=500:100:5000
    V = rand(N, N)

    V[1,:] = boundary_X_Up
    V[:,1] = boundary_Y_R
    V[N,:] = boundary_X_Down
    V[:,N] = boundary_Y_L

    V = solver!(V,n_reps)
    append!(listconvJacobi,abs(mean(V-Space)))
end
listconvJacobi