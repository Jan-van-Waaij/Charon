using DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      GZip # for zipped files

anchorpops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
truenumindvec = [8] # [1, 2, 4, 8]
allcombinations = [(anchorpop, truenumind) for anchorpop in anchorpops for truenumind in truenumindvec] # 104 elementen. 4*104=416 threads
treshold = 9 

for (anchorpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("storeliksEmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv.gz"))
    rewrittencsvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("rewriteEmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv.gz"))
    file = GZip.open(csvfile)
    writefile = GZip.open(rewrittencsvfile, "w")
    line = readline(file)
    write(writefile, line[begin:end-1]*",larger\n")
    line = readline(file)
    while line != ""
        sp = split(line, ',')
        truenumind = parse(Int, sp[4])
        estnumind = parse(Int, sp[5])
        liktruenumind = parse(Float64, sp[treshold+truenumind])
        likestnumind = parse(Float64, sp[treshold+estnumind])
        larger = liktruenumind - likestnumind
        write(writefile, line[begin:end-1]*",$larger\n")
        line = readline(file)
    end 
    close(file)
    close(writefile)
end 

mymean(x) = mean(x[end÷2:end])

for (anchorpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("rewriteEmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv.gz"))
    df = CSV.read(csvfile, DataFrame)
    gdf = groupby(df, :chainid)
    cdf = combine(gdf, :larger => mymean => :mean)
    println(cdf)
    println()
end 

for (anchorpop, truenumind) in allcombinations[2:end]
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("rewriteEmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile*".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    select!(df, Not(:larger))
    CSV.write(csvfile, df)
    command = `gzip -f $csvfile`
    run(command)
end 


for (anchorpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("storeliksEmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv.gz"))
    rewrittencsvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("rewriteEmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv.gz"))
    cp(rewrittencsvfile, csvfile)
end 



