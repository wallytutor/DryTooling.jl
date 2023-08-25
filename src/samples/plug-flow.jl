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

selected = ["CH4", "O2", "CO2", "H2O", "N2"]
gas = dry.IdealGasMixture(data, selected)
