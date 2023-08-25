# -*- coding: utf-8 -*-
# module PlugFlow

using ModelingToolkit
using Plots
using YAML

include("../src/IdealGas.jl")
using .IdealGas

function loaddatabase(fpath)
    return YAML.load_file(realpath(fpath))
end

CANTERA_DATA = "C:\\Program Files\\Cantera\\data"
data = loaddatabase(joinpath(CANTERA_DATA, "gri30.yaml"))
species = 

selected = ["CH4", "O2", "CO2", "H2O", "N2"]
gas = IdealGasMixture(data, selected)

# N₂ = IdealGasSpecies(species, "N2")


# end # (module PlugFlow)




















































































