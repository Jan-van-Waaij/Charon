using GZip 

inputfile = "data/vcffilesMads/isoaltes_4pop_1kindi.vcf.gz"
file = GZip.open(inputfile)

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
    println(counter) # 3_612_735 before rewrite. 3_595_222  after
end 

close(file)
