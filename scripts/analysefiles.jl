using StatsBase

anchorpops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, truenumind) for anchorpop in anchorpops for truenumind in truenumindvec] # 104 elementen. 

for truenumind in truenumindvec
    for anchorpop in anchorpops
        frequencyfile = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind$(truenumind)_freq_$(anchorpop).gz"
        if !isfile(frequencyfile)
            println(frequencyfile, " does not exist.")
            continue
        end 
        freqs = parsefreqfile(frequencyfile)

        readsfile = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind$(truenumind)_basecount_IBS_ind$(truenumind).gz"
        if !isfile(readsfile)
            println(readsfile, " does not exist.")
            continue
        end 
        derivedreadsvec, coveragevec = parsereadfile(readsfile, ' ')
        #@show length(freqs) length(derivedreadsvec) length(coveragevec)
        @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
        triples = zip(freqs, derivedreadsvec, coveragevec)
    end 
end 

using Data, StatsBase
truenumind=2
anchorpop = "ACB"
frequencyfile = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind$(truenumind)_freq_$(anchorpop).gz"
if !isfile(frequencyfile)
    println(frequencyfile, " does not exist.")
end 
freqs = parsefreqfile(frequencyfile)

readsfile = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind$(truenumind)_basecount_IBS_ind$(truenumind).gz"
if !isfile(readsfile)
    println(readsfile, " does not exist.")
end 
derivedreadsvec, coveragevec = parsereadfile(readsfile, ' ')
@assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

allowedindices = 0.0 .< freqs .< 1.0
freqs = freqs[allowedindices]
derivedreadsvec = derivedreadsvec[allowedindices]
coveragevec = coveragevec[allowedindices]


triples = zip(freqs, derivedreadsvec, coveragevec)  # 5_073_550
countmap(triples) # 

#= 
I think I have to change my program to your DICE-format. Because for every site where the frequency, number of derived reads and coverage is the same, the likelihood is equal. So that means many redundant computations with my format. For example /home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind2_basecount_IBS_ind2.gz has 5_073_550 lines, while in your format there are only 70_315 lines. 

Shall I keep it in the format 

ancestral_reads	   derived_reads	 	anchor_frequency 		num_sites, 

as in DICE, or do you prefer another format? 
=#