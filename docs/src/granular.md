# DryTooling.Granular

## General porous media

Modeling of geometrical characteristics of porous beds is required for including both their thermal effect or role over chemistry in chemical reactors. A classical approach used in several commercial and open source tools is that of Gunn (1978) [^Gunn1978]. In what follows we develop the ideas that lead to an analogous model which is implemented by this structure.

To build the model we will assume a reactor of constant rectangular cross-section ``{A}_{r}={b}{w}`` and volume ``{V}_{R}={b}{w}{h}``. Its cross-section perimeter is then ``{P}_{R}=2({b}+{w})``. Inside this reactor we randomly pack cubic particles ``\beta`` of surface area ``{A}_{\beta}=6{l}_{\beta}^2`` and volume ``{V}_{\beta}={l}_{\beta}^3`` at a porosity level ``{\phi}``. Thus the total volume of solids inside the reactor is ``{V}_{S}=(1-{\phi}){V}_{R}`` and the approximate number of particles ``{N}=\frac{{V}_{S}}{{V}_{\beta}}``. Following a similar reasoning the total surface area of particles is ``{A}_{S}={N}{A}_{\beta}``. Performing all the substitutions so far one finds the following expression

```math
{A}_{S}=\frac{6(1-{\phi}){b}{w}{h}}{{l}_{\beta}}
```

Since the differential ``d{A}={P}d{l}`` holds for the surface of a body over its length ``{l}``, one can divide the above expression by the reactor length to get the perimeter of particles in a cross-section. We can further divide by the cross-section area itself and find the *perimeter density* which is a more general result, and find the expression proposed by Gunn [^Gunn1978]. This result is summarized in the next equation where the subscript of particle size was dropped for generality.

```math
{P} = \frac{6(1-{\phi})}{{l}}
```

An estimator of the number of channels per unit cross-section of reactor ``{N}`` can be related to the porosity through ``{N}\pi{R}^2={phi}``. Because the above perimeter is shared between the fluid volume and solids, it holds that ``{N}2\pi{R}=P``. Using these expressions one can solve for the porosity channels characteristic *radius* ``{R}`` as given below, which is also a result reported by Gunn [^Gunn1978].

```math
{R}=\frac{{\phi}{l}}{3(1-{\phi})}
```

```@docs
DryTooling.Granular.PackedBedPorosityDescriptor
```

## Rotary kiln models

Implements the ordinary differential equation for prediction of bed height profile in a rotary kiln as proposed by Kramers and Croockewite (1952) [^Kramers1952]. Its goal is to be used as a process support tool or to integrate more complex models requiring integration of the bed profile.

```@docs
DryTooling.Granular.SymbolicLinearKramersModel
```

Description of a rotary kiln bed geometry computed from the solution of bed height along the kiln length. The main goal of the quantities computed here is their use with heat and mass transfer models for the simulation of rotary kiln process.

```@docs
DryTooling.Granular.RotaryKilnBedSolution
DryTooling.Granular.plotlinearkramersmodel
```

## References

[^Gunn1978]: [D. J. Gunn, 1978](https://doi.org/10.1016/0017-9310(78)90080-7)

[^Kramers1952]: [Kramers et al., 1952](https://doi.org/10.1016/0009-2509(52)87019-8)
