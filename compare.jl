using Data
 
freqfileGabriel = "/home/projects/gabriel/Jan_allelefreq/simulations/split/all_freq_CEU.gz"

freqsGabriel = parsefreqfile(freqfileGabriel)

freqfileJan = "data/frequencies/freqanchor/Madsmodel_freq_anchor_CEU.txt.gz"

freqsJan = parsefreqfile(freqfileJan)



readfileGabriel = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind8_basecount_CHBS_ind8.gz"

readsGabriel = parsereadfile(readfileGabriel, ' ')

readfileJan = "data/reads/Mads_model_pop_CHBS_num_ind_8_cov_64_reads_with_errors.txt.gz"

readsJan = parsereadfile(readfileJan)