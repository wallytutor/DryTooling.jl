# DryTooling.Simulation

**NOTE:** this module is *fragile* and breaking changes are still expected. It is not until all the main solvers are migrated that it will become stable. This is necessary for ensuring compatibility with all models.

## Iterative solver

The core of the iterative time-stepping solver is `step!`. This function that is described below works according to the following solution logic:

1. The `fouter!` update function is called once. This is generall where one implements the right-hand side update of the problem before stepping.

2. An initial update with `finner!` is done. Normally this is responsible by the update of matrix coefficients that are dependent on solution state.

3. If relaxation `α <= 0.0`, then the problem is treated as linear.

4. Otherwise a maximum of `M` iterations are repeated, where `fsolve!` is used to solve the problem (ofter an under-relaxation step) and is also expected to keep track of residuals.

5. Problem coefficients are updated with `finner!` if not converged.

```@docs
DryTooling.Simulation.step!
```

The outer iteration for advancing between steps is carried out by `advance!`.

```@docs
DryTooling.Simulation.advance!
```

Any model willing to implement its solution through the methods provided in this module is expected to explicity import and override the behaviour of the following methods for its own type:

```@docs
DryTooling.Simulation.fouter!
DryTooling.Simulation.finner!
DryTooling.Simulation.fsolve!
DryTooling.Simulation.timepoints
```

## Linear algebra

```@docs
DryTooling.Simulation.TridiagonalProblem
DryTooling.Simulation.solve!(::DryTooling.Simulation.TridiagonalProblem)
DryTooling.Simulation.change
```

## Residuals tracking

```@docs
DryTooling.Simulation.TimeSteppingSimulationResiduals
DryTooling.Simulation.finaliterationdata
DryTooling.Simulation.addresidual!
DryTooling.Simulation.plotsimulationresiduals
```

## Examples

Please go to the module samples [page](samples.md).
