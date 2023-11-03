# for loop in bash
for name in $(find data -name "*.vcf.gz"); do for i in 1 2 3 4 5 8 16 32 64 128 256 512; do julia  scripts/simreads.jl $name $i 20; done; done;

# get error code of program
echo $?

for name in $(find . -name "*.txt.gz"); do julia ../../scripts/simreadsrandomallele.jl $name 20; done;

for name in $(find data -name "*.vcf.gz"); do for i in 1 2 3 4 5 8 16 32 64; do julia  scripts/randomalleles.jl $name $i; done; done;

for name in $(find data -name "*.vcf.gz"); do echo $name done;

for name in $(find . -name "*.txt.gz"); do julia ../../scripts/loglikrandomallelescomplex.jl $name; done;

for name in $(find . -name "*.txt.gz"); do julia ../../scripts/loglikrandomallelescomplex.jl $name; done;

for name in $(find . -name "*.vcf.gz"); do for i in 128 256 512; do julia  ../../scripts/randomalleles.jl $name $i; done; done;

# in data/vcffiles folder
for name in $(find . -name "*.vcf.gz"); do julia ../../scripts/calcfreqanchor.jl $name; done;

for name in $(find . -name "*.txt"); do dice -2p -o "../diceoutput/${name}_output.txt" $name; done;

for name in $(find . -name "*.txt"); do echo "../diceoutput/${name}_output.txt" $name; done;
dicefiles="find *.txt"
withoutextension="${dicefiles%.txt}"