# -*- coding: utf-8 -*-
import DryTooling as dry
using ModelingToolkit
using Plots
using YAML

function loaddatabase(fpath)
    return YAML.load_file(realpath(fpath))
end

CANTERA_DATA = "C:\\Program Files\\Cantera\\data"
data = loaddatabase(joinpath(CANTERA_DATA, "gri30.yaml"))

# Create gas phase only once.
selected = ["CH4", "O2", "CO2", "H2O", "N2"]
mix = dry.IdealGasMixture(data, selected)

# Create as many simulation objects as required.
gas = dry.IdealGasSolution(mix)

@parameters P
@variables T Y[1:gas.mix.nspecies]

gas.P = P
gas.T = T
gas.Y = Y


# ρ = dry.densitymass(gas)