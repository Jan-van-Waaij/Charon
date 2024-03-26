using Plots

xs = Vector{Float64}(undef, 81)
ys = Vector{Float64}(undef, 81)
zs = Vector{Float64}(undef, 81)
vals = 0.1:0.1:0.9
let 
    count = 1
    for x in vals 
        for y in vals 
            xs[count] = x
            ys[count] = y 
            zs[count] = x*(1-x)*y*(1-y)
        end
    end 
end  

p = surface(xs, ys, zs)
savefig(p, "testplot.png")
