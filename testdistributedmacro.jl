using Distributed

addprocs(10)


@everywhere function waitandsee(filenumber)
    sleep(10)
    println(pwd())
    println(filenumber)
    file = "testtextfile$(filenumber).txt"
    openedfile = open(file, "w")
    println(openedfile, filenumber)
    close(openedfile)
    return 3 
end 

@time @sync @distributed for number in 1:8
    waitandsee(number)
end 
println("next")

