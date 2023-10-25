# DryTooling.Granular theory

```@contents
Pages = ["theory.md"]
Depth = 3
```

## Geometrical properties of granular media

Modeling of geometrical characteristics of porous beds is required for including both their thermal effect or role over chemistry in chemical reactors. A classical approach used in several commercial and open source tools is that of Gunn [Gunn1978](@cite). In what follows we develop the ideas that lead to an analogous model which is implemented by this structure.

To build the model we will assume a reactor of constant rectangular cross-section ``{A}_{r}={b}{w}`` and volume ``{V}_{R}={b}{w}{h}``. Its cross-section perimeter is then ``{P}_{R}=2({b}+{w})``. Inside this reactor we randomly pack cubic particles ``\beta`` of surface area ``{A}_{\beta}=6{l}_{\beta}^2`` and volume ``{V}_{\beta}={l}_{\beta}^3`` at a porosity level ``{\phi}``. Thus the total volume of solids inside the reactor is ``{V}_{S}=(1-{\phi}){V}_{R}`` and the approximate number of particles ``{N}=\frac{{V}_{S}}{{V}_{\beta}}``. Following a similar reasoning the total surface area of particles is ``{A}_{S}={N}{A}_{\beta}``. Performing all the substitutions so far one finds the following expression

```math
{A}_{S}=\frac{6(1-{\phi}){b}{w}{h}}{{l}_{\beta}}
```

Since the differential ``d{A}={P}d{l}`` holds for the surface of a body over its length ``{l}``, one can divide the above expression by the reactor length to get the perimeter of particles in a cross-section. We can further divide by the cross-section area itself and find the *perimeter density* which is a more general result, and find the expression proposed by Gunn [Gunn1978](@cite). This result is summarized in the next equation where the subscript of particle size was dropped for generality.

```math
{P} = \frac{6(1-{\phi})}{{l}}
```

An estimator of the number of channels per unit cross-section of reactor ``{N}`` can be related to the porosity through ``{N}\pi{R}^2={\phi}``. Because the above perimeter is shared between the fluid volume and solids, it holds that ``{N}2\pi{R}=P``. Using these expressions one can solve for the porosity channels characteristic *radius* ``{R}`` as given below, which is also a result reported by Gunn [Gunn1978](@cite).

```math
{R}=\frac{{\phi}{l}}{3(1-{\phi})}
```

## Powder bed profile in rotary drums

in a rotary kiln as proposed by Kramers and Croockewite (1952) [Kramers1952](@cite). Its goal is to be used as a process support tool or to integrate more complex models requiring integration of the bed profile. In its classical statement, the bed height profile ``h(z)`` can be evaluated from *volume* of flowing material conservation through the following equations. Coordinate ``z=0`` represents the discharge position where initial condition must be applied. This is given by the dam height, if any, or particle size.

```math
\begin{aligned}
\dfrac{dh}{dz} &= C₁ f\left(\frac{h}{R}\right) - C₂\\[6pt]
C₁             &= \frac{3}{4}\dfrac{Φ\tan{γ}}{π R^3 ω}\\[6pt]
C₂             &= \dfrac{\tan{β}}{\cos{γ}}\\[6pt]
f(r)           &= (2r - r²)^{-\frac{3}{2}}
\end{aligned}
```
