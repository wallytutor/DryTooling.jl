### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ b7781285-0b98-439f-85b4-d1b9cf72919a
begin
    @info "Loading packages"
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using CairoMakie
    using Trapz
    using PlutoUI
    import DryTooling as dry
end

# ╔═╡ 1d74deee-9948-443c-a04a-7b1f81f53bd1
include("parameters.jl")

# ╔═╡ b307aa70-641c-11ee-27cb-b151d31fffc9
md"""
# Diffusion Calculations

$(TableOfContents())

This notebook aims at validating the solution of energy transport in solids for coupling with gas phase in counter-flow model. This includes the solution of 1-D spherical particles time-dependent solution, which are used in deterministic or stochastic approaches of mean particle temperature determination, and the axial dispersion model for approximating vertical transport of energy in solids.
"""

# ╔═╡ 98022d94-143c-400e-85b9-e83fd3e9b689
md"""
## Cylinder temperature model

Heat equation formulated with temperature as dependent variable is stated as:

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=\nabla\cdotp{}(k\nabla{}T)
```

For computing the heating dynamics in a cylinder, using the definition of divergence in cylindrical coordinates and using the gradient expansion over the radius we have

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=
\frac{1}{r}\frac{\partial}{\partial{}r}
\left(rk\frac{\partial{}T}{\partial{}r}\right)
```

To proceed with the finite volume discretization we perform the integration of both sides of the equation over the relevant variables. The order of integration is chosen according to the nature of the derivative term, as discussed by Patankar (1980). Care must be taken in the definition of the space integration, which is non-trivial in cylindrical coordinates systems and must be carried over the differential volume ``dV``.

```math
\int_{V}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}dtdV=
\int_{0}^{\tau}\int_{V}
\frac{1}{r}\frac{\partial}{\partial{}r}
\left(rk\frac{\partial{}T}{\partial{}r}\right)dVdt
```

This differential volume is given by ``dV=rdr{}d\theta{}dz``. Since the problem is specified to be symmetric around cylinder center (this must include initial conditions), the azimuthal and axial components can be moved outside the time and radial integrations and lead to a common ``2\pi{}z`` factor in both sides of the equation, which cancels out.

```math
\int_{0}^{z}\int_{0}^{2\pi}d\theta{}dz=2\pi{}z
```

The integration over radial coordinate introduces the ``rdr`` factor from the differential volume and we get the final form of the equation to integrate.

```math
\int_{s}^{n}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}rdtdr=
\int_{0}^{\tau}\int_{s}^{n}
\frac{\partial}{\partial{}r}
\left(rk\frac{\partial{}T}{\partial{}r}\right)drdt
```

Effecting the inner integrations and moving out constant terms from the integrals we have

```math
\rho{}c_{p}\left(T_P^{\tau}-T_P^{0}\right)\int_{s}^{n}rdr=
\int_{0}^{\tau}
\left(rk\frac{\partial{}T}{\partial{}r}\right)\bigg\vert_{s}^{n}dt
```

Expanding the evaluation of the definite integral between control volume boundaries ``s`` and ``n`` and performing a Crank-Nicolson integration of the right-hand side one gets

```math
\begin{align}
\frac{\rho{}c_{p}}{\tau}
\left(T_P^{\tau}-T_P^{0}\right)
\left(\frac{r_n^2}{2}-\frac{r_s^2}{2}\right)
&=f\left[
r_nk_n\frac{T_N^{\tau}-T_P^{\tau}}{\delta_{P,N}}-
r_sk_s\frac{T_P^{\tau}-T_S^{\tau}}{\delta_{P,S}}
\right]\\[8pt]
&+(1-f)\left[
r_nk_n\frac{T_N^{0}-T_P^{0}}{\delta_{P,N}}-
r_sk_s\frac{T_P^{0}-T_S^{0}}{\delta_{P,S}}
\right]
\end{align}
```

Some coefficients appearing in the above equations are now grouped. Notice that for thermal conductivity ``k`` which is a function of temperature, the corresponding time-step temperature must be used for its evaluation. For ``\beta_{j}`` the lower case ``j`` represents the evaluation at the interface with control volume ``J``, what is a very specific notation.

```math
\begin{align}
\alpha_{P}  & = \frac{\rho{}c_{p}}{2\tau}\left(r_n^2-r_s^2\right)\\[8pt]
\beta_{j}   & = \frac{r_jk_j}{\delta_{P,J}}
\end{align}
```

For conciseness we make ``g=(1-f)`` and simplify the expression with the new coefficients as

```math
-f\beta_{s}T_S+
(\alpha_{P}+f\beta_{n}+f\beta_{s})T_P
-f\beta_{n}T_N
=
g\beta_{s}T_S^{0}+
(\alpha_{P}-g\beta_{n}-g\beta_{s})T_P^{0}+
g\beta_{n}T_N^{0}
```

### Implicit implementation

For the fully implicity time-stepping scheme ``f=1`` the expression reduces to

```math
-\beta_{s}T_S+
(\alpha_{P}+\beta_{n}+\beta_{s})T_P
-\beta_{n}T_N
=
\alpha_{P}T_P^{0}
```

where the following coefficients are identified

```math
\begin{align}
a_{S} & = -\beta_{s}\\[8pt]
a_{N} & = -\beta_{n}\\[8pt]
a_{P} & = \alpha_{P}+\beta_{n}+\beta_{s}
\end{align}
```

and the standard format FVM discretization is reached

```math
a_ST_S + a_PT_P + a_NT_N = \alpha_{P}T_P^{0}
```

A condition for symmetry is that no flux traverses the center of the cylinder at ``r=0``. That implies that *south* derivatives in discretizes form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
a_1T_P + a_NT_N = \alpha_{P}T_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```

Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
a_ST_S + a_RT_P = \alpha_{P}T_P^{0}+UT_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+U+\beta_{s}
```

It must be noted here that ``U=Rh``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.
"""

# ╔═╡ bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
md"""
## Sphere temperature model

In the case of spherical coordinates we start with a modification in divergence operator as follows

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)
```

The expression is again integrated over time and the differential volume ``dV``.

```math
\int_{V}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}dtdV=
\int_{0}^{\tau}\int_{V}
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)dVdt
```

This differential volume is given by ``dV=r^2\sin\phi{}dr{}d\theta{}d\phi``. Since the problem is specified to be symmetric around sphere center (this must include initial conditions), the polar and azimuthal components can be moved outside the time and radial integrations and lead to a common ``4\pi`` factor in both sides of the equation, which cancels out.

```math
\int_{0}^{\pi}\int_{0}^{2\pi}\sin\phi{}d\theta{}d\phi=4\pi
```

The integration over radial coordinate introduces the ``r^2dr`` factor from the differential volume and we get the final form of the equation to integrate.

```math
\int_{s}^{n}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}r^2dtdr=
\int_{0}^{\tau}\int_{s}^{n}
\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)drdt
```

After effecting the inner integrations and moving out constant terms from the integrals and expanding the evaluation of the definite integral between control volume boundaries ``s`` and ``n`` and performing a Crank-Nicolson integration of the right-hand side one gets

```math
\begin{align}
\frac{\rho{}c_{p}}{\tau}
\left(T_P^{\tau}-T_P^{0}\right)
\left(\frac{r_n^3}{3}-\frac{r_s^3}{3}\right)
&=f\left[
r_n^2k_n\frac{T_N^{\tau}-T_P^{\tau}}{\delta_{P,N}}-
r_s^2k_s\frac{T_P^{\tau}-T_S^{\tau}}{\delta_{P,S}}
\right]\\[8pt]
&+(1-f)\left[
r_n^2k_n\frac{T_N^{0}-T_P^{0}}{\delta_{P,N}}-
r_s^2k_s\frac{T_P^{0}-T_S^{0}}{\delta_{P,S}}
\right]
\end{align}
```

Some coefficients appearing in the above equations are now grouped. Notice that for thermal conductivity ``k`` which is a function of temperature, the corresponding time-step temperature must be used for its evaluation. For ``\beta_{j}`` the lower case ``j`` represents the evaluation at the interface with control volume ``J``, what is a very specific notation.

```math
\begin{align}
\alpha_{P}  & = \frac{\rho{}c_{p}}{3\tau}\left(r_n^3-r_s^3\right)\\[8pt]
\beta_{j}   & = \frac{r_j^2k_j}{\delta_{P,J}}
\end{align}
```

For conciseness we make ``g=(1-f)`` and simplify the expression with the new coefficients as

```math
-f\beta_{s}T_S+
(\alpha_{P}+f\beta_{n}+f\beta_{s})T_P
-f\beta_{n}T_N
=
g\beta_{s}T_S^{0}+
(\alpha_{P}-g\beta_{n}-g\beta_{s})T_P^{0}+
g\beta_{n}T_N^{0}
```

### Implicit implementation

For the fully implicity time-stepping scheme ``f=1`` the expression reduces to

```math
-\beta_{s}T_S+
(\alpha_{P}+\beta_{n}+\beta_{s})T_P
-\beta_{n}T_N
=
\alpha_{P}T_P^{0}
```

where the following coefficients are identified

```math
\begin{align}
a_{S} & = -\beta_{s}\\[8pt]
a_{N} & = -\beta_{n}\\[8pt]
a_{P} & = \alpha_{P}+\beta_{n}+\beta_{s}
\end{align}
```

and the standard format FVM discretization is reached

```math
a_ST_S + a_PT_P + a_NT_N = \alpha_{P}T_P^{0}
```

A condition for symmetry is that no flux traverses the center of the sphere at ``r=0``. That implies that *south* derivatives in discretizes form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
a_1T_P + a_NT_N = \alpha_{P}T_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```

Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
a_ST_S + a_RT_P = \alpha_{P}T_P^{0}+UT_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+U+\beta_{s}
```

It must be noted here that ``U=R^2h``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.
"""

# ╔═╡ 26ea8a50-8fb6-4d7c-a623-92911555a249
function balance(m::dry.Cylinder1DTemperatureModel, ρ, c, T₀)
    t = dry.timeaxis(m)
    r = m.grid.r
    T = m.problem.x
    Q = m.mem[].Q
    R = last(r)

    # Initial enthalpy.
    q0 = π * R^2 * ρ * c * T₀

    # Time integration => ∫₀ᵗ Q dt
    qt = trapz(t, Q)

    # Space integration => ∫₀ᴿ ρcTrdr * ∫dθ = 2πρc*∫₀ᴿTrdr
    qr = 2π * ρ * c * trapz(r, r .* T) - q0

    # Multiply by 2D to get same order of magnitude as the sphere.
    qt * 2R, qr * 2R
end

# ╔═╡ 4fefd939-d204-494e-848d-535d0f3259a3
function balance(m::dry.Sphere1DTemperatureModel, ρ, c, T₀)
    t = dry.timeaxis(m)
    r = m.grid.r
    T = m.problem.x
    Q = m.mem[].Q
    R = last(r)

    # Initial enthalpy.
    q0 = (4/3) * π * R^3 * ρ * c * T₀

    # Time integration => ∫₀ᵗ Q dt
    qt = trapz(t, Q)

    # Space integration => ∫₀ᴿ ρcTr²dr * ∫sin(ϕ)dϕdθ = 4πρc*∫₀ᴿTr²dr
    qr = 4π * ρ * c * trapz(r, r.^2 .* T) - q0

    qt, qr
end

# ╔═╡ 2ed55310-c24d-4c64-93d2-b07a852d642c
md"""
### Crank-Nicolson generalization
"""

# ╔═╡ 484ad8a3-acfe-4eda-8b79-ad2c14d6d327
md"""
## Sphere enthalpy model

Heat equation for a constant density phase using enthalpy as dependent variable is stated as:

```math
\rho{}\frac{\partial{}h}{\partial{}t}=\nabla\cdotp{}(k\nabla{}T)
```

For computing the heating dynamics in a sphere, using the definition of divergence in spherical coordinates and using the gradient expansion over the radius we have

```math
\rho{}\frac{\partial{}h}{\partial{}t}=
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)
```

This is now integrated over the differential volume ``dV`` as described in previous sections and for conciseness we skip that discussion. The integration over radial coordinate introduces the ``r^2dr`` factor from the differential volume and we get the final form of the equation to integrate.

```math
\int_{s}^{n}\int_{0}^{\tau}
\rho{}\frac{\partial{}h}{\partial{}t}r^2dtdr=
\int_{0}^{\tau}\int_{s}^{n}
\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)drdt
```

After effecting the inner integrations and moving out constant terms from the integrals and expanding the evaluation of the definite integral between control volume boundaries ``s`` and ``n`` and performing a Crank-Nicolson integration of the right-hand side one gets

```math
\begin{align}
\frac{\rho{}}{\tau}
\left(h_P^{\tau}-h_P^{0}\right)
\left(\frac{r_n^3}{3}-\frac{r_s^3}{3}\right)
&=f\left[
r_n^2k_n\frac{T_N^{\tau}-T_P^{\tau}}{\delta_{P,N}}-
r_s^2k_s\frac{T_P^{\tau}-T_S^{\tau}}{\delta_{P,S}}
\right]\\[8pt]
&+(1-f)\left[
r_n^2k_n\frac{T_N^{0}-T_P^{0}}{\delta_{P,N}}-
r_s^2k_s\frac{T_P^{0}-T_S^{0}}{\delta_{P,S}}
\right]
\end{align}
```

Some coefficients appearing in the above equations are now grouped. Notice that for thermal conductivity ``k`` which is a function of temperature, the corresponding time-step temperature must be used for its evaluation. For ``\beta_{j}`` the lower case ``j`` represents the evaluation at the interface with control volume ``J``, what is a very specific notation.

```math
\begin{align}
\alpha_{P}  & = \frac{\rho{}}{3\tau}\left(r_n^3-r_s^3\right)\\[8pt]
\beta_{j}   & = \frac{r_j^2k_j}{\delta_{P,J}}
\end{align}
```

For conciseness we make ``g=(1-f)`` and simplify the expression with the new coefficients as

```math
\begin{align}
\alpha_{P}h_P^{\tau}-\alpha_{P}h_P^{0}

&=f\beta_{n}T_N^{\tau}-f(\beta_{n}+\beta_{s})T_P^{\tau}-f\beta_{s}T_S^{\tau}
\\[8pt]
&+g\beta_{n}T_N^{0}-g(\beta_{n}+\beta_{s})T_P^{0}-g\beta_{s}T_S^{0}
\end{align}
```
"""

# ╔═╡ 66ec78ad-2154-4a53-835b-2c3a15fcdf8f
md"""
### Implicit implementation

For the fully implicity time-stepping scheme ``f=1`` and making ``\gamma_{j}=\alpha_{P}^{-1}\beta_{j}`` one gets

```math
h_P^{\tau}=h_P^{0}+\gamma_{n}T_N^{\tau}-(\gamma_{n}+\gamma_{s})T_P^{\tau}-\gamma_{s}T_S^{\tau}
```



This is no longer a linear problem and thus cannot be solved directly. We need now an strategy for solving this coupled system of nonlinear equations.
"""

# ╔═╡ 86341ac4-d32e-463e-9d19-e91ade707d06


# ╔═╡ 81312dba-9bc5-4c9c-aad1-7eac30954cfb
md"""



A condition for symmetry is that no flux traverses the center of the cylinder at ``r=0``. That implies that *south* derivatives in discretizes form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
a_1T_P + a_NT_N = \alpha_{P}T_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```

Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
a_ST_S + a_RT_P = \alpha_{P}T_P^{0}+UT_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+U+\beta_{s}
```

It must be noted here that ``U=Rh``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.
"""

# ╔═╡ 20a61a1b-0e71-4b8f-9687-d7e836a4831d
md"""
## Advection-diffusion plug-flow
"""

# ╔═╡ 1552e5db-4a33-47ca-a66f-fe4fafb40945
md"""
## Tools
"""

# ╔═╡ 883a0bef-9743-4dfd-8e63-aec333e6a5d5
"Display temperature evolution kymograph."
function temperaturekymograph(; model, cmap = :afmhot, clims = (300, 1500))
    t = dry.timeaxis(model)
    r = 100model.grid.r
    T = model.mem[].T
    R = r[end]

    fig, ax, hm = heatmap(t, r, T, colormap = cmap)

    cb = Colorbar(fig[:, end+1], hm)
    cb.limits = clims

    ax.title  = "Temperature kymograph"
    ax.xlabel = "Time [s]"
    ax.ylabel = "Radial position [cm]"

    ax.yticks = 0.0:1.0:R
    ylims!(ax, (0.0, R))

    return fig
end

# ╔═╡ 09ac1638-83a4-4b86-96ba-12d2ddd00ba6
"Plot temperature profile over sphere radius."
function plotspheretemperature(r, T, B)
    R = r[end]
    xticks = 0.0:1.0:100R
    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)
    scatter!(ax, 100R, B, color = :red)
    lines!(ax, 100r, T)
    ax.xlabel = "Radial coordinate [cm]"
    ax.ylabel = "Temperature [K]"
    ax.xticks = xticks
    xlims!(ax, (-0.2, 100R+0.2))
    return fig
end

# ╔═╡ 68383996-f87b-42b6-94e6-c8086c503d7c
ctm, figctm = let
    @info "Simulating with `Cylinder1DTemperatureModel`"
    R = 0.5blocksize
    tol = 1.0e-10

    model = dry.Cylinder1DTemperatureModel(;
        grid = dry.CylinderGrid1DEquispaced(R, 100),
        κ    = (_) -> 2.0,
        ρ    = 3000.0,
        c    = 900.0,
        h    = 20.0,
        B    = 1500.0
    )

    @time residuals = dry.solve(model;
        T     = 300.0,
        τ     = 2.0,
        t     = 2400.0,
        iters = 20,
        relax = 0.001,
        tol   = tol,
    )

    fig1 = dry.plotsimulationresiduals(residuals; ε = tol)
    fig2 = plotspheretemperature(model.grid.r, model.problem.x, model.B)
    fig3 = temperaturekymograph(; model = model)
    model, (fig1, fig2, fig3)
end;

# ╔═╡ ba35642b-f5db-4d6d-b737-4e17fdf5ecb5
figctm[1]

# ╔═╡ 297eaf10-4aa9-4c1c-90ed-f4a2b57e0e9d
figctm[2]

# ╔═╡ 26d2ae75-0298-4656-af11-a5d63d6935a8
figctm[3]

# ╔═╡ 2a3a4451-8cf8-43f9-90fb-5dd036d8895a
let
    ρ  = 3000.0
    c  = 900.0
    T₀ = 300.0
    qt, qr = balance(ctm, ρ, c, T₀)
    qt, qr, qt / qr
end

# ╔═╡ ce1817d3-5dea-4237-8be4-4770c3e08796
stm, figstm = let
    @info "Simulating with `Sphere1DTemperatureModel`"
    R = 0.5blocksize
    tol = 1.0e-10

    model = dry.Sphere1DTemperatureModel(;
        grid = dry.SphereGrid1DEquispaced(R, 100),
        κ    = (_) -> 2.0,
        ρ    = 3000.0,
        c    = 900.0,
        h    = 20.0,
        B    = 1500.0
    )

    @time residuals = dry.solve(model;
        T     = 300.0,
        τ     = 2.0,
        t     = 2400.0,
        iters = 20,
        relax = 0.001,
        tol   = tol,
    )

    fig1 = dry.plotsimulationresiduals(residuals; ε = tol)
    fig2 = plotspheretemperature(model.grid.r, model.problem.x, model.B)
    fig3 = temperaturekymograph(; model = model)
    model, (fig1, fig2, fig3)
end;

# ╔═╡ 12a1a5e9-57fe-40e9-a039-e96dc8e4a9bf
figstm[1]

# ╔═╡ ddb18f9f-58c2-4512-b8a6-43d032fb629e
figstm[2]

# ╔═╡ ef9aaecd-4c40-4e37-8d89-179a3cd67b63
figstm[3]

# ╔═╡ d5154b1b-1c99-4a3a-93cc-74fede3830c9
let
    ρ  = 3000.0
    c  = 900.0
    T₀ = 300.0
    qt, qr = balance(stm, ρ, c, T₀)
    qt, qr, qt / qr
end

# ╔═╡ Cell order:
# ╟─b307aa70-641c-11ee-27cb-b151d31fffc9
# ╟─98022d94-143c-400e-85b9-e83fd3e9b689
# ╟─68383996-f87b-42b6-94e6-c8086c503d7c
# ╟─ba35642b-f5db-4d6d-b737-4e17fdf5ecb5
# ╟─297eaf10-4aa9-4c1c-90ed-f4a2b57e0e9d
# ╟─26d2ae75-0298-4656-af11-a5d63d6935a8
# ╟─bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
# ╟─ce1817d3-5dea-4237-8be4-4770c3e08796
# ╟─12a1a5e9-57fe-40e9-a039-e96dc8e4a9bf
# ╟─ddb18f9f-58c2-4512-b8a6-43d032fb629e
# ╟─ef9aaecd-4c40-4e37-8d89-179a3cd67b63
# ╟─26ea8a50-8fb6-4d7c-a623-92911555a249
# ╟─4fefd939-d204-494e-848d-535d0f3259a3
# ╟─2a3a4451-8cf8-43f9-90fb-5dd036d8895a
# ╟─d5154b1b-1c99-4a3a-93cc-74fede3830c9
# ╟─2ed55310-c24d-4c64-93d2-b07a852d642c
# ╟─484ad8a3-acfe-4eda-8b79-ad2c14d6d327
# ╠═66ec78ad-2154-4a53-835b-2c3a15fcdf8f
# ╠═86341ac4-d32e-463e-9d19-e91ade707d06
# ╠═81312dba-9bc5-4c9c-aad1-7eac30954cfb
# ╟─20a61a1b-0e71-4b8f-9687-d7e836a4831d
# ╟─1552e5db-4a33-47ca-a66f-fe4fafb40945
# ╟─883a0bef-9743-4dfd-8e63-aec333e6a5d5
# ╟─09ac1638-83a4-4b86-96ba-12d2ddd00ba6
# ╟─b7781285-0b98-439f-85b4-d1b9cf72919a
# ╟─1d74deee-9948-443c-a04a-7b1f81f53bd1
