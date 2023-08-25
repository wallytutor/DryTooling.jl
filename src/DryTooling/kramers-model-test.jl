# -*- coding: utf-8 -*-
import DryTooling as dry
using Statistics

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
