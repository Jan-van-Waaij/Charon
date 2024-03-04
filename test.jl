println(Base.JLOptions().opt_level)

println(Threads.nthreads())

println(Threads.nthreads() == Sys.CPU_THREADS)

using Distributed
println(nprocs())
