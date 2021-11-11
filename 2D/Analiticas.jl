using Plots

L=2
function φ(x,y;V=10)
    L=2
    s=0
    for n=1:2:100000
        s=s.+map(z  -> sin(n *π *z/L) ,x).*map( z -> exp(-n *π *z/L)./n,y)
    end
    4s*V/ π
end


Space = zeros(100,100)
x=0:2/(99):2;y=x


length(x)
for i ∈ 1:100  , j ∈ 1:100
    Space[i,j]=φ(x[i],y[j])
end
plot(x,y,Space,linetype=:contourf,title = "Analitica")
