using GZip 

inputfile = "data/vcffiles_demes/multiple_demes_multipledemes_12.vcf.gz"


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
    println(counter) # 61767 before rewrite. 61632 after
end 

close(file)
