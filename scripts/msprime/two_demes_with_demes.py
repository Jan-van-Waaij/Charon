import sys # for command line arguments. 
import msprime 
import demes
import subprocess # to call gzip from within Python. 

def two_demes(demesfile, demesfiledescription, random_seed=None):
    graph = demes.load(demesfile)
    print('1')
    demography = msprime.Demography.from_demes(graph)
    print('2')
    ts = msprime.sim_ancestry(samples={"ancient":5_000, "anchor":5_000}, 
    demography=demography, random_seed=random_seed, recombination_rate=1e-8, 
    sequence_length=10e7)
    print('3')
    mts = msprime.sim_mutations(ts, rate=2e-8, random_seed=random_seed)
    print('4')
    filenamewithoutextension = "two_demes_with_demes_" + \
    demesfiledescription + '_' + str(random_seed)
    print('5')
    mts.dump( filenamewithoutextension+".tree")
    print('6')
    vcffilename=filenamewithoutextension+".vcf"
    print('7')
    with open(vcffilename, "w") as vcf_file:
        mts.write_vcf(vcf_file)
    print('8')
    subprocess.call("/home/ctools/bin/bgzip -f "+vcffilename, shell=True)
    print('9')

if len(sys.argv) == 4:
    two_demes( sys.argv[1], sys.argv[2], int(sys.argv[3]))
elif len(sys.argv) == 3:
    two_demes(sys.argv[1], sys.argv[2])
else: 
    println("Something went wrong, use two or three arguments: the first argument should be the demes file, the second one a description, and the third is optionally the seed.")