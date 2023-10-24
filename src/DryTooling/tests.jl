# -*- coding: utf-8 -*-
using DryTooling

@testset "Basic functionalities" begin
    @test closestpowerofx(10) ≈ 10.0
    @test closestpowerofx(11) ≈ 20.0
    @test closestpowerofx(11; roundf = floor) ≈ 10.0
    @test closestpowerofx(11; x = 5, roundf = floor) ≈ 10.0

    v = collect(1:4)
    @test head(v) == [1; 2; 3]
    @test tail(v) == [2; 3; 4]
    @test body(v) == [2; 3]

    @test heaviside(-1) == 0
    @test heaviside(-1.0) == 0.0
    @test heaviside(0.0) == 0.5
    @test heaviside(1.0) == 1.0
    @test interval(10; a = 0, b = 10) == 0.5

    # TODO makestepwise1d
end
