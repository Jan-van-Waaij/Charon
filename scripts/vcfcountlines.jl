using GZip 

inputfile = "data/vcffiles_demes/two_demes_with_demes_tauC=0.25tauA=0.25numgen=1000effpop=200012.vcf.gz"

file = GZip.open(inputfile)
println("hi")

let 
    line = readline(file)
    counter = 0
    while line != ""
        counter += 1
        line = readline(file)
        if counter % 100_000 == 0
            println(counter)
        end 
    end 
    println(counter) # 28_287_219 before rewrite. 28_260_464 after
end 

close(file)
