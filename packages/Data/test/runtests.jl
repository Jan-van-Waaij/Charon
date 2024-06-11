using Test, Data, Distributions, GZip

@testset "Data" begin
    @testset "DataError" begin
        @test_throws Data.DataError throw(Data.DataError())
        @test_throws Data.DataError throw(Data.DataError(1,2))
        @test_throws Data.DataError throw(Data.DataError("something wrong"))
        @test_throws Data.DataError throw(Data.DataError("something ", "wrong"))
    end

    @testset "FileError" begin
        file = "bla.txt"
        @test_throws Data.FileError throw(Data.FileError())
        @test_throws Data.FileError throw(Data.FileError("something wrong with $file."))
        @test_throws Data.FileError throw(Data.FileError("something wrong"))
        @test_throws Data.FileError throw(Data.FileError("something ", "wrong"))
    end

    # skiplines(openedfile::IO, numberoflinestobeskipped::Int)
    @testset "skiplines" begin
        GZip.open("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/1to10.txt.gz") do file 
            skiplines(file, 0)
            line = readline(file)
            number = parse(Int, line)
            @test number == 1
        end
        GZip.open("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/1to10.txt.gz") do file
            skiplines(file, 7)
            line = readline(file)
            number = parse(Int, line)
            @test number == 8
        end 
        GZip.open("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/1to10.txt.gz") do file
            @test_throws DomainError skiplines(file, -1)
            @test_throws DomainError skiplines(file, -10)
        end 
    end

    @testset "isderived" begin
        @test Data.isderived('1')
        @test !Data.isderived('0')
        @test !Data.isderived('|')
        @test_throws Data.DataError Data.isderived('a')
    end

    # countderivedalleles(line::Vector{SubString{String}}, population::AbstractVector{Int})
    @testset "countderivedalleles" begin
        sp = split("0|0\t1|0\t0|0\t1|0\t1|0\t0|1\t0|0\t1|1") # 8 genotypes, 
        @test Data.countderivedalleles(sp, [1, 3]) == 0
        @test Data.countderivedalleles(sp, 1:8) == 6
        # @test_throws exception expr
        sp = split("0|0\t1|0\t0|0\t1|0\t1|0\t0|1\t0|0\t1|1\tbla") # wrong format
        @test_throws Data.DataError Data.countderivedalleles(sp, [9])
    end

    @testset "addhashtagstoeachline" begin 
        @test addhashtagstoeachline("") == ""
        text = """
        bla 
        bla bla 
        bla """
        textwithhashtags = addhashtagstoeachline(text)
        sp = split(textwithhashtags, '\n')
        sp = sp[begin:end-1] # last item is "", because
        # the last character of textwithhashtags is \n. 
        @test length(sp) == 3
        for line in sp 
            @test line[1] == '#'
        end 
    end 

    @testset "openfile" begin
        @test_throws Data.FileError Data.openfile("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/1to10.txt.gz", ".vcf.gz")
        openedfile = Data.openfile("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/1to10.txt.gz", ".txt.gz")
        @test close(openedfile) == 0 # hack to see if file was open. 
    end

    @testset "addcomments" begin
        f = GZip.open("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/tmp.txt.gz", "w")
        comments = """
        bla 
        bla 
        bla """
        Data.addcomments(f, comments)
        close(f)
        f = GZip.open("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/tmp.txt.gz")
        line = readline(f)
        while line != ""
            @test line[1] == '#'
            line = readline(f)
        end 
        close(f)
        rm("/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/tmp.txt.gz")
    end

    @testset "openwritefileandwritecomments" begin
        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/tmp.txt.gz"
        comments = """
        bla 
        bla 
        bla """
        openedfile = Data.openwritefileandwritecomments(file, false, comments)
        close(openedfile)
        openedfile = GZip.open(file)
        line = readline(openedfile)
        while line != ""
            @test line[1] == '#'
            line = readline(openedfile)
        end 
        close(openedfile)
        isfile(file) || error("make sure the file exist for testing")
        @test_throws Data.FileError Data.openwritefileandwritecomments(file, false, comments)
        openedfile = Data.openwritefileandwritecomments(file, true, comments)
        close(openedfile)
        openedfile = GZip.open(file)
        line = readline(openedfile)
        while line != ""
            @test line[1] == '#'
            line = readline(openedfile)
        end 
        close(openedfile)
        rm(file)
    end

    @testset "reads" begin
        @test reads(10, 10, 5)[2] == 5
        @test reads(10, 10, 10)[2] == 10
        numexp = 10_000 
        derreads = Vector{Int}(undef, numexp)
        coverages = Vector{Int}(undef, numexp)
        for index in eachindex(derreads, coverages)
            derreads[index], coverages[index] = reads(5,10,100)
        end 
        @test unique(coverages) == [100]
        @test abs(mean(derreads) - 25.0) < 1.0
        for index in eachindex(derreads, coverages)
            derreads[index], coverages[index] = reads(5,10,Poisson(100))
        end 
        @test abs(mean(coverages) - 100.0) < 1.0
        @test abs(mean(derreads) - 25.0) < 1.0

        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        openedfile = GZip.open(file)
        individualsinthesediment = sample(1:400, 10, replace = false)
        coverage = 20
        allderivedreads, allcoverages = reads(openedfile, individualsinthesediment, 20)
        close(openedfile)
        @test unique(allcoverages) == [20]
        @test abs(mean(allderivedreads) - coverage/2) < 1.0

        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        individualsinthesediment = sample(1:400, 10, replace = false)
        coverage = 20
        allderivedreads, allcoverages = reads(file, individualsinthesediment, 20)
        @test unique(allcoverages) == [20]
        @test abs(mean(allderivedreads) - coverage/2) < 1.0

        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        openedfile = GZip.open(file)
        writefile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/temp.vcf.gz"
        openedwritefile = GZip.open(writefile, "w")    
        individualsinthesediment = sample(1:400, 10, replace = false)
        coverage = 20
        reads(openedfile, openedwritefile, individualsinthesediment, 20)
        close(openedfile)
        close(openedwritefile)
        openedwritefile = GZip.open(writefile)
        allderivedreadsfile, allcoveragesfile = parsereadfile(openedwritefile)
        close(openedwritefile)
        rm(writefile)
        @test unique(allcoveragesfile) == [20]
        @test abs(mean(allderivedreadsfile) - coverage/2) < 1.0

        # from here with Poisson
        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        openedfile = GZip.open(file)
        individualsinthesediment = sample(1:400, 10, replace = false)
        coverage = 20
        allderivedreads, allcoverages = reads(openedfile, individualsinthesediment, Poisson(20))
        close(openedfile)
        @test abs(mean(allcoverages) - 20) < 1.0
        @test abs(mean(allderivedreads) - coverage/2) < 1.0

        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        individualsinthesediment = sample(1:400, 10, replace = false)
        coverage = 20
        allderivedreads, allcoverages = reads(file, individualsinthesediment, Poisson(20))
        @test abs(mean(allcoveragesfile) -20) < 1.0
        @test abs(mean(allderivedreadsfile) - coverage/2) < 1.0

        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        openedfile = GZip.open(file)
        writefile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/temp.vcf.gz"
        openedwritefile = GZip.open(writefile, "w")    
        individualsinthesediment = sample(1:400, 10, replace = false)
        coverage = 20
        reads(openedfile, openedwritefile, individualsinthesediment, Poisson(20))
        close(openedfile)
        close(openedwritefile)
        openedwritefile = GZip.open(writefile)
        allderivedreadsfile, allcoveragesfile = parsereadfile(openedwritefile)
        close(openedwritefile)
        rm(writefile)
        @test abs(mean(allcoveragesfile) -20) < 1.0
        @test abs(mean(allderivedreadsfile) - coverage/2) < 1.0
    end

    @testset "calc_freq" begin
        line = "0|0\t1|1\t0|1\t1|0\t1|1\t1|1\t0|0"#0 2 1 1 2 2 0
        @test Data.calc_freq(line, 1:7) == 8/14
        @test Data.calc_freq(line, 5:7) == 2/3
        @test Data.calc_freq(line, 2:4) == 2/3
    end

    @testset "freq" begin
        # with opened files. 
        vcffile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        openedfile = GZip.open(vcffile)
        f = freq(openedfile, 1:400)
        close(openedfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        freqfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/temp.txt.gz"
        openvcffile = GZip.open(vcffile)
        openfreqfile = GZip.open(freqfile, "w")
        freq(openvcffile, openfreqfile, 1:400)
        close(openvcffile)
        close(openfreqfile)
        openfreqfile = GZip.open(freqfile)
        f = parsefreqfile(openfreqfile)
        close(openfreqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000

        openvcffile = GZip.open(vcffile)
        openfreqfile = GZip.open(freqfile, "w")
        freq(openvcffile, openfreqfile, 1:400, "comments\nmorecomments\n")
        close(openvcffile)
        close(openfreqfile)
        openfreqfile = GZip.open(freqfile)
        f = parsefreqfile(openfreqfile)
        close(openfreqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000

        openvcffile = GZip.open(vcffile)
        openfreqfile = GZip.open(freqfile, "w")
        freq(openvcffile, openfreqfile, 1:400, "comments\nmorecomments")
        close(openvcffile)
        close(openfreqfile)
        openfreqfile = GZip.open(freqfile)
        f = parsefreqfile(openfreqfile)
        close(openfreqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000

        # with strings. 
        vcffile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        f = freq(vcffile, 1:400)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        freqfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/temp.txt.gz"
        rm(freqfile, force=true) # remove freqfile if it exists. 
        freq(vcffile, freqfile, 1:400)
        f = parsefreqfile(freqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        rm(freqfile, force=true) # remove freqfile if it exists.
        freq(vcffile, freqfile, 1:400, "comments\nmorecomments\n")
        f = parsefreqfile(freqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        rm(freqfile, force=true) # remove freqfile if it exists.
        freq(vcffile, freqfile, 1:400, "comments\nmorecomments")
        f = parsefreqfile(freqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000

        # with strings overwrite file 
        vcffile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/fiftypercentderived.vcf.gz"
        f = freq(vcffile, 1:400)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        freqfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/temp.txt.gz"
        freq(vcffile, freqfile, 1:400, true)
        f = parsefreqfile(freqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        freq(vcffile, freqfile, 1:400, true, "comments\nmorecomments\n")
        f = parsefreqfile(freqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000
        freq(vcffile, freqfile, 1:400, true, "comments\nmorecomments")
        f = parsefreqfile(freqfile)
        @test unique(f) == [0.5]
        @test length(f) == 10_000

        isfile(freqfile) || error("freqfile does not exist. Needed for proper testing.")
        @test_throws Data.FileError freq(vcffile, freqfile, 1:400)
        @test_throws Data.FileError freq(vcffile, freqfile, 1:400, "comments\nmorecomments")
        @test_throws Data.FileError freq(vcffile, freqfile, 1:400, false)
        @test_throws Data.FileError freq(vcffile, freqfile, 1:400, false, "comments\nmorecomments")
    end

    @testset "parsefreqfile" begin
        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/freqfile.txt.gz" # file with 11 lines n-th lines
        # 0.0
        # 0.1 
        # ...
        # 1.0 
        isfile(file) || error("$file does not exist.")
        allfreqs = Data.parsefreqfile(file)
        @test allfreqs == 0.0:0.1:1.0
        openedfile = GZip.open(file)
        allfreqs = Data.parsefreqfile(openedfile)
        @test allfreqs == 0.0:0.1:1.0
        close(openedfile)
    end

    @testset "parsereadfile" begin
        file = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/packages/Data/test/reads.txt.gz" # file with 11 lines n-th lines
        # 0 10
        # 1 10 
        # ...
        # 5 10
        # 6 20
        # 7 20
        # ...
        # 10 20
        isfile(file) || error("$file does not exist.")
        derivedreads, coverages = Data.parsereadfile(file)
        @test derivedreads == 0:10
        @test coverages == vcat(repeat([10],6), repeat([20],5))
        openedfile = GZip.open(file)
        allfreqs = Data.parsereadfile(openedfile)
        @test derivedreads == 0:10
        @test coverages == vcat(repeat([10],6), repeat([20],5))
        close(openedfile)
    end

    @testset "subsample" begin 
        @test subsample(3, 5, 5) == 3, 5
        @test subsample(5, 5, 3) == 3
        @test subsample(0, 5, 3) == 0
        for oldcoverage in 2:10
            for newcoverage in 1:oldcoverage
                for derivedreads in 0:oldcoverage 
                    @test subsample(derivedreads, oldcoverage, newcoverage)[1] ≤ derivedreads
                end 
            end 
        end 
        results = Vector{Int}(undef, 1000_000)
        for rate in 0.1:0.1:0.9
            for oldcoverage in 1:100
                for newcoverage in 1:oldcoverage
                    for index in eachindex(results)
                        derivedreads = rand(Binomial(oldcoverage, rate))
                        results[index] = subsample(derivedreads, oldcoverage, newcoverage)
                    end
                    newrate = mean(results)/newcoverage
                    @test abs(rate-newrate) < 0.01
                end
            end
        end
    end # subsample testset. 
end 