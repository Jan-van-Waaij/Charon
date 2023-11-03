# rundice.jl 
# Author: Jan van Waaij
# Created: unknown, but before 16 February 2023.
# Last update: 16 February 2023.
# 
# Purpose:
# Run dice. 
#
# Instalation:
# Please read the Julia package install instructions in installinstructions.txt.
#
# Usage:
# julia pathto/rundice.jl 
# 
# Output:
# Ouput files are saved in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/diceoutputs
# as filename_dice_output.txt, for each filename_dicefile.txt in 
# /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/dicefiles. 

path = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/dicefiles"
dicefiles = readdir(path)
dicefileending = "_dicefile.txt"
isdicefile = (x -> endswith(x, dicefileending))
filter!(isdicefile, dicefiles)
dicefilepaths = joinpath.(path, dicefiles)

dicefileendinglength = length(dicefileending)
dicefileoutputs = [file[1:end-dicefileendinglength]*"_dice_output.txt" for file in dicefiles]
diceoutputspath = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/diceoutputs"
dicefileoutputsfullpaths = joinpath.(diceoutputspath, dicefileoutputs)
commands = similar(dicefileoutputsfullpaths, Cmd)
for index in eachindex(commands, dicefileoutputsfullpaths, dicefilepaths)
    input = dicefilepaths[index]
    output = dicefileoutputsfullpaths[index]
    commands[index] = `/home/ctools/dice/src/dice -2p -o $output $input`
end

run.(commands)
