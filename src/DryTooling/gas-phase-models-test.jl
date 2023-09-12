# -*- coding: utf-8 -*-
import DryTooling as dry
using ModelingToolkit

# @testset "Gas Phase Models" begin
#     @test begin
#        selected = ["CH4", "O2", "CO2", "H2O", "N2"]
#        gas = dry.IdealGasMixture(data, selected)
#        M = gas.molecularmasses
#        species = gas.species[end]
#
#        Tnum = 1000.0
#        @parameters Tsym
#        @variables Tvar
#
#        Ynum = ones(gas.nspecies) / gas.nspecies
#        @parameters Ysym[1:gas.nspecies]
#
#        gas.Y[1:end] = Ynum
#
#        dry.specificheatmass(species, Tnum)
#        dry.specificheatmass(species, Tsym)
#        dry.specificheatmass(species, Tvar)
#
#        dry.enthalpymass(species, Tnum)
#        dry.enthalpymass(species, Tsym)
#        dry.enthalpymass(species, Tvar)
#
#        dry.specificheatmole(species, Tnum)
#        dry.specificheatmole(species, Tsym)
#        dry.specificheatmole(species, Tvar)
#
#        dry.enthalpymole(species, Tnum)
#        dry.enthalpymole(species, Tsym)
#        dry.enthalpymole(species, Tvar)
#
#        dry.meanmolecularmass(M, Ynum)
#        dry.meanmolecularmass(M, Ysym)
#
#        dry.massfraction2molefraction(M, Ynum)
#        dry.massfraction2molefraction(M, Ysym)
#
#        dry.molefraction2massfraction(M, Ynum)
#        dry.molefraction2massfraction(M, Ysym)
#
#        dry.meanmolecularmass(gas)
#        dry.densitymass(gas)
#        sum(dry.massfractions(gas)) ≈ 1.0
#        sum(dry.molefractions(gas)) ≈ 1.0
#        dry.specificheatmass(gas)
#     end
# end