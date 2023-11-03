using GZip 

inputfile = "data/vcffiles_simulated_by_Jan/two_demes_time_in_generations_1000_2000_500_1.vcf.gz"
outputfile = "data/vcffiles_simulated_by_Jan/two_demes_time_in_generations_1000_2000_500_1_rewrite.vcf.gz"
input = GZip.open(inputfile)
output = GZip.open(outputfile, "w")

for i in 1:6
    line = readline(input)
    print(output, line)
end 

line = readline(input)
while line != ""
    sp = split(line)
    if unique(sp[10:end]) ⊆ ["0|0", "0|1", "1|0", "1|1"]
        print(output, line)
    end 
    global line = readline(input)
end 

close(input)
close(output)

rm(inputfile)
sleep(10)
cp(outputfile, inputfile)
sleep(10)
rm(outputfile)
