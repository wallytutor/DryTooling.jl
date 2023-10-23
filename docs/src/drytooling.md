# DryTooling Core

## Physical constants

```@docs
DryTooling.GAS_CONSTANT
DryTooling.ZERO_CELSIUS
DryTooling.ONE_ATM
DryTooling.STEFAN_BOLTZMANN
```

## Haskell-like array slicing

```@docs
DryTooling.head
DryTooling.tail
DryTooling.body
```

## Handling of discontinuous functions

```@docs
DryTooling.heaviside
DryTooling.interval
DryTooling.makestepwise1d
```

## Rounding numbers and automatic axes

```@docs
DryTooling.closestpowerofx
DryTooling.axesunitscaler
```

## Abstract types

### Problem solving and physical models

```@docs
DryTooling.AbstractMatrixProblem
DryTooling.AbstractIterativeSolver
DryTooling.AbstractSolutionStorage
DryTooling.AbstractPhysicalModel
```

### Transport, thermodynamics, and kinetics

```@docs
DryTooling.AbstractTransportModel
DryTooling.AbstractSolidTransport
DryTooling.AbstractGasThermo
DryTooling.AbstractSolidThermo
DryTooling.AbstractSolidMaterial
DryTooling.AbstractMixtureSubstance
DryTooling.AbstractMixturePhase
DryTooling.AbstractKineticsMechanism
```

### Finite volume method and relatives

```@docs
DryTooling.AbstractDiffusionModel1D
DryTooling.AbstractGrid1D
```

## Other constants

```@docs
DryTooling.STABLE_ELEMENTS_TABLE
DryTooling.TRANSPORT_MODELS
```
