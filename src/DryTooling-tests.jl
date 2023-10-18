# -*- coding: utf-8 -*-
import DryTooling as dry
using ModelingToolkit
using Statistics

@testset "Basic functionalities" begin
    @test dry.closestpowerofx(10) ≈ 10.0
    @test dry.closestpowerofx(11) ≈ 20.0
    @test dry.closestpowerofx(11; roundf = floor) ≈ 10.0
    @test dry.closestpowerofx(11; x = 5, roundf = floor) ≈ 10.0

    v = collect(1:4)
    @test dry.head(v) == [1; 2; 3]
    @test dry.tail(v) == [2; 3; 4]
    @test dry.body(v) == [2; 3]

    @test dry.heaviside(-1) == 0
    @test dry.heaviside(-1.0) == 0.0
    @test dry.heaviside(0.0) == 0.5
    @test dry.heaviside(1.0) == 1.0
    @test dry.interval(10; a = 0, b = 10) == 0.5

    # TODO makestepwise1d
end

@testset "Grid generation" begin
end

@testset "Simulation Residuals" begin
end

@testset "Gas Phase Models" begin
    # @test begin
    #    selected = ["CH4", "O2", "CO2", "H2O", "N2"]
    #    gas = dry.IdealGasMixture(data, selected)
    #    M = gas.molecularmasses
    #    species = gas.species[end]

    #    Tnum = 1000.0
    #    @parameters Tsym
    #    @variables Tvar

    #    Ynum = ones(gas.nspecies) / gas.nspecies
    #    @parameters Ysym[1:gas.nspecies]

    #    gas.Y[1:end] = Ynum

    #    dry.specificheatmass(species, Tnum)
    #    dry.specificheatmass(species, Tsym)
    #    dry.specificheatmass(species, Tvar)

    #    dry.enthalpymass(species, Tnum)
    #    dry.enthalpymass(species, Tsym)
    #    dry.enthalpymass(species, Tvar)

    #    dry.specificheatmole(species, Tnum)
    #    dry.specificheatmole(species, Tsym)
    #    dry.specificheatmole(species, Tvar)

    #    dry.enthalpymole(species, Tnum)
    #    dry.enthalpymole(species, Tsym)
    #    dry.enthalpymole(species, Tvar)

    #    dry.meanmolecularmass(Ynum, M)
    #    dry.meanmolecularmass(Ysym, M)

    #    dry.massfraction2molefraction(Ynum, M)
    #    dry.massfraction2molefraction(Ysym, M)

    #    dry.molefraction2massfraction(Ynum, M)
    #    dry.molefraction2massfraction(Ysym, M)

    #    dry.meanmolecularmass(gas)
    #    dry.densitymass(gas)
    #    sum(dry.massfractions(gas)) ≈ 1.0
    #    sum(dry.molefractions(gas)) ≈ 1.0
    #    dry.specificheatmass(gas)
    # end
end

@testset "Shomate material" begin
    # Webbook NIST ID=C14808607&Type=JANAF&Table=on
    silica = dry.MaterialShomate(
        a_lo = [-6.076591, 251.6755, -324.7964, 168.5604,
                0.002548, -917.6893, -27.96962, -910.8568],
        a_hi = [58.7534, 10.27925, -0.131384, 0.02521,
                0.025601, -929.3292, 105.8092, -910.8568],
        T_ch = 847.0
    )
    T = [298.0, 300.0, 400.0, 847.0, 900.0, 1900.0]
    c = [44.57, 44.77, 53.43, 67.42, 67.95, 77.99]
    @test sum(abs2, silica.cₚ.(T) - c) < 0.0001
end

@testset "Kramers Model" begin
    R = 1.0e+00
    Φ = 1.0e-02
    z = collect(0.0:0.1:10.0)

    # Force a very thick *analytical* bed.
    h = R * ones(size(z))

    # Analytical *half-filled* bed area.
    Aₐ = π * R^2 / 2

    # Analytical *half-filled* bed volume.
    Vₐ = Aₐ * z[end]

    # Create bed object.
    bed = dry.RotaryKilnBedSolution(z, h, R, Φ)

    @test mean(bed.θ) ≈ π
    @test mean(bed.l) ≈ 2R
    @test mean(bed.A) ≈ Aₐ
    @test mean(bed.η) ≈ 0.5
    @test bed.ηₘ      ≈ 50.0
    @test bed.V       ≈ Vₐ
    @test bed.τ       ≈ Vₐ / Φ
end
