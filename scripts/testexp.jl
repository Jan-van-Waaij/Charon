using Likelihood, PrecompileTools

Qs = [Matrix(makeQ(n)) for n in 1:30]
eigenQs = EigenExpansion.(Qs)
Qꜜs = [Matrix(makeQꜜ(n)) for n in 1:30]
eigenQꜜs = EigenExpansion.(Qꜜs)

@compile_workload begin
    calcmatrix(Qs[1], Qꜜs[1], 0.3, 0.2)
    calcmatrix(eigenQs[1], eigenQꜜs[1], 0.3, 0.2)
    calcmatrix(Ref(eigenQs, 1)[], Ref(eigenQꜜs, 1)[], 0.3, 0.2)
    calcmatrix(view(eigenQs, 1)[], view(eigenQꜜs, 1)[], 0.3, 0.2)
end

println("Time old version")
@time for n in 1:30
        for τC in 0.01:0.01:1.00
            for τA in 0.01:0.01:1.00
                calcmatrix(Qs[n], Qꜜs[n], τC, τA)
            end
        end 
    end

println("Time new version with getindex")
@time for n in 1:30
        for τC in 0.01:0.01:1.00
            for τA in 0.01:0.01:1.00
                calcmatrix(eigenQs[n], eigenQꜜs[n], τC, τA)
            end
        end 
    end

println("Time new version with Ref")
@time for n in 1:30
        for τC in 0.01:0.01:1.00
            for τA in 0.01:0.01:1.00
                calcmatrix(Ref(eigenQs, n)[], Ref(eigenQꜜs, n)[], τC, τA)
            end
        end 
    end

println("Time new version with view()[]")
@time for n in 1:30
        for τC in 0.01:0.01:1.00
            for τA in 0.01:0.01:1.00
                calcmatrix(view(eigenQs, n)[], view(eigenQꜜs, n)[], τC, τA)
            end
        end 
    end

# De nieuwe versies zijn bijna allemaal ongeveer even snel.
# view(A, index)[] ietsje pietje sneller in twee keer testen
# maar die heeft wel meer allocations. 
