using Distributed 

node = [("compute01", 20)]

myworkers = addprocs(node)

@everywhere begin
    using Distributed
    myid 
    sleep 
end 

println("I start!")
@sync @distributed for _ in 1:20 
    sleep(60*10)
end  

#=
println("I start!")
@time begin 
    try 
        @sync begin  
            @distributed for _ in 1:20 
                sleep(60*10)
            end
            println("everyone is sleeping")
        end  
    catch e
        #kill.(myworkers)
        interrupt()
        println("I killed all workers!")
    end 
end  
=#