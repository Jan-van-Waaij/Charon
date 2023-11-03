using PosteriorSampling, Data, CSV, DataFrames, Plots

using Data, # For parsing frequency and read files. 
      GZip, # Read and write zip files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates # for the now function 

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "mcmcnolocalmaximacoverage100.csv.gz")

openedcsvfile = GZip.open(csvfile, "a") # the a of append. 

df = DataFrame(NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    trueτC = Float64[],
                    trueτA = Float64[],
                    truenumind = Int[],
                    coverage = Int[],
                    nsample = Int[],
                    τCsample = Float64[],
                    τAsample = Float64[],
                    accepted = Bool[],
                    chainid = Int[]
                    )

CSV.write(openedcsvfile, df)

filewithoutextension = "two_demes_11000_9000_70000_0"

frequencyfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqanchor/two_demes_11000_9000_70000_0_freqanchor.txt.gz"

y = parsefreqfile(frequencyfile)

R = 100

readsfile = string("/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads/two_demes_11000_9000_70000_0_4_", R, "_reads.txt.gz")



d = parsereadfile(readsfile)[1]

sp = split(filewithoutextension, '_')
NeAncientPop = parse(Int, sp[3])
NeAnchorPop = parse(Int, sp[4])
jointtimeyears = parse(Int, sp[5])
generationtime = 25
trueτA = jointtimeyears / (2*NeAncientPop*generationtime)
trueτC = jointtimeyears / (2*NeAnchorPop*generationtime)


prioronn = DiscreteUniform(1, 8)
prioronτCτA = Product([Uniform(), Uniform()])
nsteps = 100_000

results = MHsampler(4, nsteps, prioronn, prioronτCτA, R, d, y) 
for i in eachindex(results) 
    nsample, τCsample, τAsample, acceptedvec = results[i]
    local df = DataFrame(NeAnchorPop = repeat([NeAnchorPop], nsteps),
                NeAncientPop = repeat([NeAncientPop], nsteps),
                jointtimeyears = repeat([jointtimeyears], nsteps),
                trueτC = repeat([trueτC], nsteps),
                trueτA = repeat([trueτA], nsteps),
                truenumind = repeat([4], nsteps),
                coverage = repeat([R], nsteps),
                nsample = nsample,
                τCsample = τCsample,
                τAsample = τAsample,
                accepted = acceptedvec,
                chainid = repeat([Int(i)], nsteps))
    CSV.write(openedcsvfile, df, append=true)
end # i in results.

close(openedcsvfile)
