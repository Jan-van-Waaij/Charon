# Created 13 - 10 - 2023
# Set of Julia instructions to analyse the exact posterior.
# Meant to use interactively. 

using CSV, DataFrames

CSV.write("data/fullposterior/exactposteriorsediment1Xunknownn.csv", df)

df = CSV.read("data/fullposterior/exactposteriorsediment1Xunknownn.csv", DataFrame)

df[:, Not(:localmaxima, :localmaximavalues)]

df.localmaxima = eval.(Meta.parse.(df.localmaxima))

df.localmaximavalues = eval.(Meta.parse.(df.localmaximavalues))

show(df[:, [:localmaxima]], truncate=200)

