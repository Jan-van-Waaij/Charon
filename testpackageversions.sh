for csvversion in {1..10}; do 
    csvversionstring="0.$csvversion"
    for distributionsversion in {1..25}; do
        distributionsversionstring="0.$distributionsversion" 
        for statsbaseversion in {1..34}; do
            statsbaseversionstring="0.$statsbaseversion"
            rm -f Testversions/Manifest.toml
            rm -f Testversions/Project.toml

            julia testpackageversions.jl $csvversionstring $distributionsversionstring $statsbaseversionstring
        done 
    done
done 
