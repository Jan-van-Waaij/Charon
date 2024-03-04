using GZip 

inputfile = "/home/projects/MAAG/msprime_deme/janvan/isoaltes_4pop_1kindi.vcf.gz"
outputfile = "data/vcffilesMads/isoaltes_4pop_1kindi.vcf.gz"
input = GZip.open(inputfile)
output = GZip.open(outputfile, "w")

for i in 1:6 # copy the comments without change. 
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