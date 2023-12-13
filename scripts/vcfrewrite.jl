using GZip 

inputfile = "data/vcffiles_demes/multiple_demes_multipledemes_12.vcf.gz"
outputfile = "data/vcffiles_demes/multiple_demes_multipledemes_12_rewrite.vcf.gz"
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
cp(outputfile, inputfile)
rm(outputfile)
