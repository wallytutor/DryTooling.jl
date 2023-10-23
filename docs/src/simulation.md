# DryTooling.Simulation

**NOTE:** this module is *fragile* and breaking changes are still expected.
It is not until all the main solvers are migrated that it will become stable.
This is necessary for ensuring compatibility with all models.

## Residuals tracking

```@docs
DryTooling.Simulation.TimeSteppingSimulationResiduals
DryTooling.Simulation.finaliterationdata
DryTooling.Simulation.addresidual!
DryTooling.Simulation.plotsimulationresiduals
```

## Tutorial: residuals tracking in a solver

The residuals tracking functionalities of module `DryTooling.Simulation`` are not often imported by the end-user (except for its plotting utility function). In this tutorial we illustrate the logic of using a residual tracker in a new solver.

```julia
using DryTooling.Simulation

N = 2      # Number of variables.
M = 5      # Maximum inner steps.
steps = 10 # Time-advancement steps.

# Create a TimeSteppingSimulationResiduals object with the number of variables
# to track, how many inner iterations per step are expected, and the
# number of steps.
#
# IMPORTANT: If the total number of iterations is exceeded, it is up
# to the user to allocate more memory, the tracker will not manage it!
r = TimeSteppingSimulationResiduals(N, M, steps)

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
# constructor that accepts a `TimeSteppingSimulationResiduals` object. This
# handles the post-processing.
s = TimeSteppingSimulationResiduals(r)

# The new object is ready for visualization. Check the documentation
# of the following function for more details. It provides a raw figure
# and handles for modifying it for proper display (uncomment below).
# fig = plotsimulationresiduals(s; showinner = true)[1]
```
