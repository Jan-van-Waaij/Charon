using Distributed

nodes = vcat(["node0$x" for x in 1:8], ["node$x" for x in 11:14], ["compute01"])

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--threads=4"]

processes = Int[]
for node in nodes 
    append!(processes, addprocs([node], dir=dir, exename=exename, exeflags=exeflags))
end 

@everywhere using Hwloc 

futures = [@spawnat process num_physical_cores() for process in processes] 

cores = fetch.(futures)

@everywhere Base.Threads 

futures5 = [@spawnat process Threads.nthreads() for process in processes]
num_threads = fetch.(futures5)

for (node, numcores) in zip(nodes, cores)
    println(node, " : ", numcores)
end 

futures2 = [@spawnat process Base.active_project() for process in processes]

active_projects = fetch.(futures2)

futures3 = [@spawnat process num_virtual_cores() for process in processes]
virtual_cores = fetch.(futures3)

for (node, numcores) in zip(nodes, virtual_cores)
    println(node, " : ", numcores)
end 

maxthreadspernode = Int.(round.(virtual_cores .* 0.75, RoundDown))
maxprocessespernode = maxthreadspernode .÷ 4 
sum(maxprocessespernode)

for (node, maxprocs) in zip(nodes, maxprocessespernode)
    println(node, " : ", maxprocs)
end 

#=
node01 : 4 leeg 
node02 : 3
node03 : 4
node04 : 4
node05 : 4
node06 : 4
node07 : 12
node08 : 3
node11 : 7
node12 : 3
node13 : 3
node14 : 3
compute01 : 24
=#

@everywhere using Conda, PyCall 
Conda.add("psutil")

@everywhere py"""
import psutil
def get_cpu():
    return psutil.cpu_percent()
"""

@everywhere function getcpuusage()
    x = py"get_cpu()"
    return x 
end 

futures4 = [@spawnat process getcpuusage() for process in processes] 
fetch.(futures4)

# behalve node05, zijn is voor alle nodes het aantal virtual cores het dubbele van het aantal physical nodes. 
# Bij node05 is het aantal virtual nodes hetzelfde als het aantal physical nodes. 

