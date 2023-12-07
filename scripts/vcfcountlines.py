import gzip

inputfile = "data/vcffiles_demes/two_demes_with_demes_tauC=0.25tauA=0.25numgen=1000effpop=200012.vcf.gz"

f = gzip.open(inputfile)

line = f.readline()

counter = 0
while (line != b""):
    counter += 1 
    line = f.readline()
    if (counter % 100_000 == 0):
        print(counter)

print(counter)