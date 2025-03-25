file = "compatibleversionsonlynumbers.txt"

df =CSV.read(file, DataFrame, header=false, names=["CSV", "Distributions", "StatsBase"])

CSVrange = []
Distributionsrange = []
StatsBaserange = []

function rowisindataframe(df, row)
    sizedf= size(df)
    sizedf[2] == length(row) || throw(DimensionMismatch("row should have the same length as the number of columns in the dataframe."))
    for n in 1:sizedf[1]
        if Vector(df[n, :]) == row 
            return true 
        end 
    end 
    return false 
end 

function checkcompetability(df, CSVrange, Distributionsrange, StatsBaserange)
    for CSVversion in CSVrange 
        for Distributionsversion in Distributionsrange 
            for StatsBaseversion in StatsBaserange
                if !rowisindataframe(df, [CSVversion, Distributionsversion, StatsBaseversion])
                    return false 
                end 
            end 
        end 
    end 
    return true 
end 

#=
Wat is een goede range van versies 
1. Nieuwste versies zitten erin. 
2. 
=#

CSVrange = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
Distributionsrange = [0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25 ]
StatsBaserange = [   0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34]
checkcompetability(df, CSVrange, Distributionsrange, StatsBaserange)

println(length(CSVrange)*length(Distributionsrange)*length(StatsBaserange)) # 60. 