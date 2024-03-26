using Distributed, Hwloc

nodes = [
    "node01",
    "node02",
    "node03",
    "node04",
    "node05",
    "node06",
    "node07",
    "node08",
    "node11",
    "node12",
    "node13",
    "node14",
    "compute01"
]

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir"] 

processes = Int[]
for node in nodes 
    append!(processes, addprocs([node], dir=dir, exename=exename, exeflags=exeflags))
end 

@everywhere function cpuusage(number::Int)
    outputfile = "top-output$number.txt"
    command = pipeline(`top -b -n 1 -H`; stdout = outputfile)
    run(command)

    totalpercentageused = 0.0
    file = open(outputfile)
    line = readline(file)
    while line != ""
        line = readline(file)
    end 
    line = readline(file) # now you have column names. 
    @assert split(line)[9] == "%CPU"
    line = readline(file)
    while line != ""
        sp = split(line)
        totalpercentageused += parse(Float64, sp[9])
        line = readline(file)
    end 
    close(file)
    rm(outputfile)
    return totalpercentageused / 100.0
end 

@everywhere function top()
    command = `top`
    run(command)
end



futures = [@spawnat process cpuusage(process) for process in processes]
usages = fetch.(futures)
@everywhere using Hwloc
futures = [@spawnat process num_virtual_cores() for process in processes]
num_cores = fetch.(futures)


for (node, usage, num_core) in zip(nodes, usages, num_cores)
    println(node, ' ', usage, " cores in use of ", num_core)
end 

