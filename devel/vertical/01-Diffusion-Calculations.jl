### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ b7781285-0b98-439f-85b4-d1b9cf72919a
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using CairoMakie
    using Trapz
    using PlutoUI
    include("parameters.jl")
end

# ╔═╡ b307aa70-641c-11ee-27cb-b151d31fffc9
md"""
# Diffusion Calculations

$(TableOfContents())

This notebook aims at validating the solution of energy transport in solids for coupling with gas phase in counter-flow model. This includes the solution of 1-D spheric particles time-dependent solution, which are used in deterministic or stochastic approaches of mean particle temperature determination, and the axial dispersion model for approximating vertical transport of energy in solids.
"""

# ╔═╡ bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
md"""
## Sphere temperature model

Heat equation formulated with enthalpy as dependent variable is stated as:

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=\nabla\cdotp{}(k\nabla{}T)
```

For computing the heating dynamics in a sphere, using the definition of divergent in spheric coordinates and using the gradient expansion over the radius we have

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)
```

To proceed with the finite volume discretization we perform the integration of both sides of the equation over the relevant variables. The order of integration is chosen according to the nature of the derivative term, as discussed by Patankar (1980). Care must be taken in the definition of the space integration, which is non-trivial in spheric coordinates systems and must be carried over the differential volume ``dV``.

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
"""

# ╔═╡ 9b5b60ef-a62b-400e-8726-ade2d01b4fe6
md"""
Effecting the inner integrations and moving out constant terms from the integrals we have

```math
\rho{}c_{p}\left(T_P^{\tau}-T_P^{0}\right)\int_{s}^{n}r^2dr=
\int_{0}^{\tau}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)\bigg\vert_{s}^{n}dt
```

Expanding the evaluation of the definite integral between control volume boundaries ``s`` and ``n`` and performing a Crank-Nicolson integration of the right-hand side one gets

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
"""

# ╔═╡ e76fc7f9-fb1e-477f-876a-dc7fc5f63824
md"""
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
"""

# ╔═╡ e1bf1aa4-3f18-4a47-9c04-abedf66b78fb
md"""
### Implicit implementation

For the fully implicity time-stepping scheme ``f=1`` and the expression reduces to

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
"""

# ╔═╡ abc61e2c-0d39-4e32-9c15-dcd273b193aa
md"""
A condition for symmetry is that no flux traverses the center of the sphere at ``r=0``. That implies that *south* derivatives in discretizes form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
a_1T_P + a_NT_N = \alpha_{P}T_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```
"""

# ╔═╡ 585f2546-260f-4f63-b3e5-747be4e419bd
md"""
Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
a_ST_S + a_RT_P = \alpha_{P}T_P^{0}+UT_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+U+\beta_{s}
```

It must be noted here that ``U=4\pi{}R^2h``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.

!!! warn
    Conservation analysis is showing an error of 4π! It is probably here!

"""

# ╔═╡ 9eae615a-c10f-40ed-a590-bf328abfe6ea
"Non-linear iteration updater for model."
function spheretemperatureinner!(model::SphereTemperatureModel)::Nothing
    # Aliases to approach mathermatical formulation.
    a_p = model.problem.A.d
    a_s = model.problem.A.dl
    a_n = model.problem.A.du

    # Wall interpolated thermal conductivities.
    kp = @. model.k(model.problem.x)
    ks = kp[1:end-1]
    kn = kp[2:end]
    κ = @. 2 * ks * kn / (ks + kn)

    # Update temperature dependency of β.
    βₖ = @. κ * model.β

    # Main, lower and upper diagonal elements.
    a_s[1:end] = -βₖ
    a_n[1:end] = -βₖ
    a_p[1:end] = model.α

    # Update main but leave boundaries off.
    @. a_p[2:end-1] += βₖ[2:end] + βₖ[1:end-1]

    # Boundary conditions on main.
    a_p[1] += βₖ[1]
    a_p[end] += βₖ[end] + model.U

    return nothing
end

# ╔═╡ ce1817d3-5dea-4237-8be4-4770c3e08796
"Time-step dependent updater for model."
function spheretemperatureouter!(
        model::SphereTemperatureModel,
        nouter::Int64,
        ts::Float64
    )::Nothing
    # Follow surface heat flux and store partial solutions.
    model.Q[nouter] = model.U * (model.T∞ - model.problem.x[end])
    model.T[nouter, 1:end] = model.problem.x[1:end]

    # Problem right-hand side.
    @. model.problem.b[1:end] = model.α * model.problem.x

    # Apply boundary condition.
    model.problem.b[end] += model.U * model.T∞

    return nothing
end

# ╔═╡ 2ed55310-c24d-4c64-93d2-b07a852d642c
md"""
### Crank-Nicolson generalization
"""

# ╔═╡ 484ad8a3-acfe-4eda-8b79-ad2c14d6d327
md"""
## Sphere enthalpy model
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
    t = 0:model.τ:model.t
    z = 100model.z
    R = z[end]

    fig, ax, hm = heatmap(t, z, model.T, colormap = cmap)

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
function plotspheretemperature(z, T, T∞)
    R = z[end]
    xticks = 0.0:1.0:100R

    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)

    scatter!(ax, 100R, T∞, color = :red)
    lines!(ax, 100z, T)

    ax.xlabel = "Radial coordinate [cm]"
    ax.ylabel = "Temperature [K]"
    ax.xticks = xticks
    xlims!(ax, (-0.2, 100R+0.2))

    return fig
end

# ╔═╡ d75feb8b-9722-436c-8c5b-92d71d6b926d
stm, figstm = let
    @info "Simulating with `SphereTemperatureModel`"
    tol = 1.0e-10

    model = SphereTemperatureModel(;
        N  = 100,
        R  = 0.5blocksize,
        ρ  = 3000.0,
        c  = 900.0,
        h  = 20.0,
        T∞ = 1500.0,
        T₀ = 300.0,
        k  = (T) -> 2.0,
        t  = 1200.0,
        M  = 2*600
    )

    @time residuals = relaxationouterloop(;
        model        = model,
        updaterouter = spheretemperatureouter!,
        updaterinner = spheretemperatureinner!,
        tend         = model.t,
        tau          = model.τ,
        iters        = 20,
        relax        = 0.00001,
        tol          = tol,
        metric       = maxrelativevariation
    )

    fig1 = plotresiduals(residuals; ε = tol)
    fig2 = plotspheretemperature(model.z, model.problem.x, model.T∞)
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
begin
    # XXX: include this conservation check in model!
    R  = 0.05
    ρ  = 3000.0
    c  = 900.0
    T₀ = 300.0
    V  = (4/3) * π * R^3
    m  = ρ * V

    t = 0:stm.τ:stm.t
    r = stm.z
    T = stm.problem.x

    # ∫₀ᵗ Q dt => this is already over the surface area!
    qn = trapz(t, 0.5 * (stm.Q[1:end-1] .+ stm.Q[2:end]))

    # ∫₀ᴿ ρcTr²dr * ∫sin(ϕ)dϕdθ = 4πρc*∫₀ᴿTr²dr => the area factor is there!
    qa = 4π * ρ * c * trapz(r, r.^2 .* T) - m * c * T₀

    # TODO why 4pi??!?!?!?!
    qn, qa, qa/ (4π * qn)
end

# ╔═╡ Cell order:
# ╟─b307aa70-641c-11ee-27cb-b151d31fffc9
# ╟─bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
# ╟─9b5b60ef-a62b-400e-8726-ade2d01b4fe6
# ╟─e76fc7f9-fb1e-477f-876a-dc7fc5f63824
# ╟─e1bf1aa4-3f18-4a47-9c04-abedf66b78fb
# ╟─abc61e2c-0d39-4e32-9c15-dcd273b193aa
# ╟─585f2546-260f-4f63-b3e5-747be4e419bd
# ╟─9eae615a-c10f-40ed-a590-bf328abfe6ea
# ╟─ce1817d3-5dea-4237-8be4-4770c3e08796
# ╠═d75feb8b-9722-436c-8c5b-92d71d6b926d
# ╟─12a1a5e9-57fe-40e9-a039-e96dc8e4a9bf
# ╟─ddb18f9f-58c2-4512-b8a6-43d032fb629e
# ╟─ef9aaecd-4c40-4e37-8d89-179a3cd67b63
# ╠═d5154b1b-1c99-4a3a-93cc-74fede3830c9
# ╟─2ed55310-c24d-4c64-93d2-b07a852d642c
# ╟─484ad8a3-acfe-4eda-8b79-ad2c14d6d327
# ╟─20a61a1b-0e71-4b8f-9687-d7e836a4831d
# ╟─1552e5db-4a33-47ca-a66f-fe4fafb40945
# ╟─883a0bef-9743-4dfd-8e63-aec333e6a5d5
# ╟─09ac1638-83a4-4b86-96ba-12d2ddd00ba6
# ╟─b7781285-0b98-439f-85b4-d1b9cf72919a
