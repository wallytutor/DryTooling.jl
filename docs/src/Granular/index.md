# DryTooling.Granular

## Documentation

### General porous media

```@docs
DryTooling.Granular.PackedBedPorosityDescriptor
```

### Rotary kiln models

The structure `SymbolicLinearKramersModel` implements the Kramers' ordinary differential equation for prediction of bed height profile in a rotary kiln. This equation is implemented under the formalism of `ModelingToolkit`.

```@docs
DryTooling.Granular.SymbolicLinearKramersModel
```

Description of a rotary kiln bed geometry computed from the solution of bed height along the kiln length. The main goal of the quantities computed here is their use with heat and mass transfer models for the simulation of rotary kiln process.

```@docs
DryTooling.Granular.RotaryKilnBedSolution
DryTooling.Granular.plotlinearkramersmodel
```

Finally a set of basic equations provided for process analysis.

```@docs
DryTooling.Granular.sullivansηₘ
DryTooling.Granular.dimlessNΦ
DryTooling.Granular.dimlessNₖ
DryTooling.Granular.perrayresidence
DryTooling.Granular.kramersnlapprox
```

## Examples

Please go to the module samples [page](samples.md).

## Theory guide

Please go to the module theory guide [page](theory.md).

## Models validation

- [Kramers' model](validation/kramers-model.md)
