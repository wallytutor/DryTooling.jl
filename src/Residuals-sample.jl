# -*- coding: utf-8 -*-
"""
DryTooling.Residuals sample
===========================

Module DryTooling.Residuals is not generally imported by the end-user except
for its plotting utility function. In this script we illustrate the logic of
using a residual tracker in a new solver, thus it is seem more as a tutorial.
"""

import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using DryTooling.Residuals

begin
    # Number of variables.
    N = 2

    # Maximum inner steps.
    M = 5

    # Time-advancement steps.
    steps = 10

    # Create a SimulationResiduals object with the number of variables
    # to track, how many inner iterations per step are expected, and the
    # number of steps.
    #
    # IMPORTANT: If the total number of iterations is exceeded, it is up
    # to the user to allocate more memory, the tracker will not manage it!
    r = SimulationResiduals(N, M, steps)

    # The following loop represents a *dummy solver*. The outer loop
    # provides the time-advancement while the inner loop handles the
    # nonlinear problem. In the inner loop we use a random number
    # generator to provide varying number of steps per outer step.
    for kouter in 1:steps
        for kinner in 1:rand(2:M)
            # Keep track of inner iterations per step.
            r.innersteps[kouter] = kinner

            # Feed residuals to the solver.
            addresidual!(r, rand(r.N))
        end
    end

    # After running a simulation we create a new object using another
    # constructor that accepts a `SimulationResiduals` object. This
    # handles the post-processing.
    s = SimulationResiduals(r)

    # The new object is ready for visualization. Check the documentation
    # of the following function for more details. It provides a raw figure
    # and handles for modifying it for proper display.
    fig = plotsimulationresiduals(s; showinner = true)[1]
end

