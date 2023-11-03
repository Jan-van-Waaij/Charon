import sys
import msprime
import subprocess

def two_demes(effectivePop1, effectivePop2, joinTimeYears, randomseed):
    # Convert the join time from years to generations.
    # The script assumes a 25-year generation time.
    generation_time = 25
    joinTime = joinTimeYears / generation_time

    # Define the demographic model.
    # Start by creating an empty demographic model.
    demography = msprime.Demography()

    # Add two populations to the demographic model with names "pop1" and "pop2",
    # and set their initial sizes.
    demography.add_population(name="pop1", initial_size=effectivePop1)
    demography.add_population(name="pop2", initial_size=effectivePop2)
    
    # Define demographic events.
    # First, at the specified join time, change the population size of "pop2"
    # to be the maximum of the initial sizes of the two populations.
    demography.add_population_parameters_change(
        time=joinTime, initial_size=max(effectivePop1, effectivePop2), growth_rate=0, population="pop2"
    )

    # Set the migration rate to 0 at the join time.
    demography.add_migration_rate_change(time=joinTime, rate=0)

    # Merge population "pop1" into "pop2" at the join time.
    demography.add_population_split(time=joinTime, derived=["pop1"], ancestral="pop2")

    # Run the coalescent simulation with the defined demography.
    # Specify the number of samples from each population, sequence length,
    # recombination rate, random seed, and ploidy.
    ts = msprime.sim_ancestry(
        samples={"pop1": 10000, "pop2": 10000},
        demography=demography,
        sequence_length=10e6,
        recombination_rate=1e-8,
        random_seed=randomseed,
        ploidy=2
    )

    # Add mutations to the simulated coalescent tree.
    # Specify the mutation rate and random seed.
    ts = msprime.sim_mutations(
        tree_sequence=ts,
        rate=2e-8,
        random_seed=randomseed
    )

    # Save the tree sequence to a file.
    # The filename includes the input parameters for easy identification.
    tree_file = f"two_demes_{effectivePop1}_{effectivePop2}_{joinTimeYears}_{randomseed}.tree"
    ts.dump(tree_file)

    # Save the simulated data as a VCF file.
    # The filename also includes the input parameters for easy identification.
    vcf_filename = f"two_demes_{effectivePop1}_{effectivePop2}_{joinTimeYears}_{randomseed}.vcf"
    with open(vcf_filename, "w") as vcf_file:
        ts.write_vcf(vcf_file, ploidy=2)

    # Compress the VCF file using bgzip.
    # This part assumes that bgzip is installed and available in /home/ctools/bin/.
    # You might need to change this depending on your system configuration.
    subprocess.call(f"/home/ctools/bin/bgzip -f {vcf_filename}", shell=True)


# Run the function with command-line arguments.
# The command-line arguments are effective population sizes of the two populations,
# the time at which they join in years, and a random seed.
two_demes(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))

