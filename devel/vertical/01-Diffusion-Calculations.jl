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

	using DryTooling: timepoints
	using DryTooling.Grids
	using DryTooling.HeatConduction
	using DryTooling.Residuals
    import DryTooling as dry

	set_theme!(
		figure_padding = 20,
		backgroundcolor = :white
	)
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
"""

# ╔═╡ bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
md"""
## Sphere temperature model
"""

# ╔═╡ 26ea8a50-8fb6-4d7c-a623-92911555a249
function balance(m::Cylinder1DTemperatureModel, ρ, c, T₀)
    t = timepoints(m)
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
function balance(m::Sphere1DTemperatureModel, ρ, c, T₀)
    t = timepoints(m)
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

### Implicit implementation

!!! warning

    Check boundary condition for inconsistency! Maybe division by α missing!

For the fully implicity time-stepping scheme ``f=1`` and making ``\gamma_{j}^{k}=\alpha_{P}^{-1}\beta_{j}^{k}`` one gets

```math
h_P^{\tau}-h_P^{0}-\gamma_{n}^{k}T_N^{\tau,k}+(\gamma_{n}^{k}+\gamma_{s}^{k})T_P^{\tau,k}-\gamma_{s}^{k}T_S^{\tau,k}=0
```

A condition for symmetry is that no flux traverses the center of the sphere at ``r=0``. That implies that *south* derivatives in discretizes form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
h_P^{\tau}-h_P^{0}-\gamma_{n}^{k}T_N^{\tau,k}+\gamma_{n}^{k}T_P^{\tau,k}=0
```

Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
h_P^{\tau}-h_P^{0}-UT_{\infty}+(U+\gamma_{s}^{k})T_P^{\tau,k}-\gamma_{s}^{k}T_S^{\tau,k}=0
```

It must be noted here that ``U=R^2h``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.

This is no longer a linear problem and thus cannot be solved directly. We need now an strategy for solving this coupled system of nonlinear equations. The iterative solution of the problem is indicated in the above equations through the introduction of superscript ``k`` indicating the iteration number. One can rework the system as

```math
\begin{align}
-\gamma_{1,2}^{k}T_2^{\tau,k}+\gamma_{1,2}^{k}T_1^{\tau,k}+h_1^{\tau}&=h_1^{0}\\
&\dots \\
-\gamma_{n}^{k}T_N^{\tau,k}+(\gamma_{n}^{k}+\gamma_{s}^{k})T_P^{\tau,k}-\gamma_{s}^{k}T_S^{\tau,k}+h_P^{\tau}&=h_P^{0}\\
&\dots \\
(U+\gamma_{K-1,K}^{k})T_K^{\tau,k}-\gamma_{K-1,K}^{k}T_{K-1}^{\tau,k}+h_K^{\tau}&=h_K^{0}+UT_{\infty}
\end{align}
```

It is clear now that for implementation purposes one can store the required coefficients in a tridiagonal matrix ``A^{k}``. Making ``\Gamma_{i}=(\gamma_{i-1,i}+\gamma_{i,i+1})`` we can identify the terms in

```math
\begin{pmatrix}
H_{1}^{k}    \\
H_{2}^{k}    \\
H_{3}^{k}    \\
\vdots   \\
H_{K-1}^{k}  \\
H_{K}^{k}    \\
\end{pmatrix}
=
\begin{pmatrix}
 \gamma_{1,2}^{k} & -\gamma_{1,2}^{k} &  0                & \dots  & 0 & 0 \\
-\gamma_{1,2}^{k} &  \Gamma_{2}^{k}   & -\gamma_{2,3}^{k} & \dots  & 0 & 0 \\
 0 & -\gamma_{2,3}^{k} &  \Gamma_{3}^{k} & -\gamma_{3,4}^{k}\ddots &  0 &  0 \\
\vdots  & \ddots & \ddots & \ddots & \ddots  & \vdots \\
 0 &  0 & 0 & -\gamma_{K-2,K-1}^{k} &  \Gamma_{K-1}^{k}   & -\gamma_{K-1,K}^{k} \\
 0      &  0     &  0     &  0     & -\gamma_{K-1,K}^{k} & U+\gamma_{K-1,K}^{k} \\
\end{pmatrix}
\begin{pmatrix}
T_{1}^{\tau,k}   \\
T_{2}^{\tau,k}   \\
T_{3}^{\tau,k}   \\
\vdots           \\
T_{K-1}^{\tau,k} \\
T_{N}^{\tau,k}   \\
\end{pmatrix}
```

Since the temperature vector ``T^{\tau,k}`` is updated every iteration, the coefficients of ``A^{k}`` must also be updated. With the intermediate vector ``H^{\tau,k}`` the nonlinear problem is rewriten as

```math
\begin{pmatrix}
H_{1}^{k}    \\
H_{2}^{k}    \\
H_{3}^{k}    \\
\vdots       \\
H_{K-1}^{k}  \\
H_{K}^{k}    \\
\end{pmatrix}
+
\begin{pmatrix}
h_{1}^{\tau}   \\
h_{2}^{\tau}   \\
h_{3}^{\tau}   \\
\vdots         \\
h_{K-1}^{\tau} \\
h_{K}^{\tau}   \\
\end{pmatrix}
=
\begin{pmatrix}
h_1^{0}                 \\
h_2^{0}                 \\
h_3^{0}                 \\
\vdots                  \\
h_{K-1}^{0}             \\
h_{K}^{0} + UT_{\infty} \\
\end{pmatrix}
```

The choice not to write the problem in this format reflects the fact that the term ``H^{\tau,k}`` on the left-hand side is updated on a iteration basis, while the vector ``b^{0}`` is computed once per time step. This last vector was called ``b^{0}`` instead of ``h^{0}`` because it also includes the boundary condition in its last element. This is useful for the conception of the inner and outer loop functions used for solution update.
"""

# ╔═╡ 795bb843-1e8a-4151-a289-532f6282ada3
md"""
The traditional approach to solve this sort of problems is to provide a *initial guess* ``T^{\tau,0}=T^{0}``.

```math
\begin{align}
h^{\tau,0}               &= b^{0}-A^{0}T^{\tau,0}\\
h(T^{\tau,1})-h^{\tau,0} &= 0\\
\Delta{}T                &= T^{\tau,1}-T^{\tau,0}\\
T^{\tau,1}               &= T^{\tau,0}+(1-\alpha)\Delta{}T\\
\varepsilon^{1}          &= \vert\Delta{}T\vert\\
&\text{repeat}\\
h^{\tau,1}               &= b^{0}-A^{1}T^{\tau,1}\\
h(T^{\tau,2})-h^{\tau,1} &= 0\\
\Delta{}T                &= T^{\tau,2}-T^{\tau,1}\\
T^{\tau,2}               &= T^{\tau,1}+(1-\alpha)\Delta{}T\\
\varepsilon^{2}          &= \vert\Delta{}T\vert\\
&\dots\\
h^{\tau,k}                 &= b^{0}-A^{k}T^{\tau,k}\\
h(T^{\tau,k+1})-h^{\tau,k} &= 0\\
\Delta{}T                  &= T^{\tau,k+1}-T^{\tau,k}\\
T^{\tau,k+1}               &= T^{\tau,k}+(1-\alpha)\Delta{}T\\
\varepsilon^{k+1}          &= \vert\Delta{}T\vert\\
\end{align}
```
"""

# ╔═╡ 294e84f7-d230-4de6-825b-b5a32068c2dd
"Thermal diffusion in a sphere represented in temperature space."
struct Sphere1DEnthalpyModel <: dry.AbstractDiffusionModel1D
    "Grid over which problem will be solved."
    grid::dry.AbstractGrid1D

    "Memory for model linear algebra problem."
    problem::dry.TridiagonalProblem

    "Constant part of model coefficient α."
    α′::Vector{Float64}

    "Constant part of model coefficient β."
    β′::Vector{Float64}

    "Thermal conductivity in terms of temperature."
    κ::Function

    "Enthalpy in terms of temperature."
    h::Function

    "Global heat transfer coefficient ``U=hR²``."
    U::Float64

    "Surface environment temperature."
    B::Float64

    "Time-step used in integration."
    τ::Base.RefValue{Float64}

    "Memory storage for solution retrieval."
    mem::Base.RefValue{dry.Temperature1DModelStorage}

    function Sphere1DEnthalpyModel(;
            grid::dry.AbstractGrid1D,
            h::Function,
            κ::Function,
            ρ::Float64,
            u::Float64,
            B::Float64
        )
        problem = dry.TridiagonalProblem(grid.N)

        rₙ = dry.tail(grid.w)
        rₛ = dry.head(grid.w)
        α′ = @. ρ * (rₙ^3 - rₛ^3) / 3.0

        rₙ = dry.tail(grid.r)
        rₛ = dry.head(grid.r)
        wⱼ = dry.body(grid.w)
        β′ = @. wⱼ^2 / (rₙ - rₛ)

        U = u * last(grid.r)^2
        τ = Ref(-Inf)
        mem = Ref(dry.Temperature1DModelStorage(0, 0))

        return new(grid, problem, α′, β′, κ, h, U, B, τ, mem)
    end
end

# ╔═╡ 0c1e97e7-5293-4504-950f-0635d62b6092
"Time-step dependent updater for model."
function sphereenthalpyouter!(
        m::Sphere1DEnthalpyModel,
        t::Float64,
        n::Int64
    )::Nothing
    # Note the factor 4π because U = r²h only and A = 4πr²!!!!
    m.mem[].Q[n] = 4π * m.U * (m.B - last(m.problem.x))
    m.mem[].T[n, 1:end] = m.problem.x

    @. m.problem.b[1:end] = m.h.(m.problem.x)
    m.problem.b[end] += m.U * m.B

    return nothing
end

# ╔═╡ 86ad23ee-f11d-45cd-8ef3-30de6afa0568
"Non-linear iteration updater for model."
function sphereenthalpyinner!(
        m::Sphere1DEnthalpyModel,
        t::Float64,
        n::Int64
    )::Nothing
    a_p = m.problem.A.d
    a_s = m.problem.A.dl
    a_n = m.problem.A.du
    T_p = m.problem.x

    κ = interfaceconductivity1D(m.κ.(T_p))
    β = κ .* m.β′
    α = m.α′./ m.τ

    a_s[1:end] = -β
    a_n[1:end] = -β
    a_p[1:end] = α

    a_p[2:end-1] += tail(β) + head(β)
    a_p[1]       += first(β)
    a_p[end]     += last(β) + m.U

    return nothing
end

# ╔═╡ f5f89e4e-8dd8-4269-8405-4e8833cfab24
function sphereenthalpysolve!()
end

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
    t = timepoints(model)
    r = 100model.grid.r
    T = model.mem[].T
    R = r[end]

    fig, ax, hm = heatmap(t, r, T, colormap = cmap)

    cb = Colorbar(fig[:, end+1], hm)
    cb.limits = clims

    ax.title  = "Temperature kymograph"
    ax.xlabel = "Time [s]"
    ax.ylabel = "Position [cm]"

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
    ε = 1.0e-10

    model = Cylinder1DTemperatureModel(;
        grid = equidistantcellsgrid1D(R, 100),
        h    = 20.0,
        B    = 1500.0,
        κ    = 2.0,
        ρ    = 3000.0,
        c    = 900.0
    )

    @time dry.solve(model;
        T     = 300.0,
        τ     = 2.0,
        t     = 2400.0,
        α     = 0.1,
        ε     = ε,
        M     = 20
    )

	residuals = SimulationResiduals(model.res[])

    fig1, ax1, p1 = plotsimulationresiduals(residuals; ε = ε)
	ax1.xticks = 0:3:12
	ax1.yticks = -11.2:0.2:-9.8
	xlims!(ax1, extrema(ax1.xticks.val))
	ylims!(ax1, extrema(ax1.yticks.val))
	
    fig2 = plotspheretemperature(model.grid.r, model.problem.x, model.B(1000000))
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
    ε = 1.0e-10

    model = Sphere1DTemperatureModel(;
        grid = equidistantcellsgrid1D(R, 100),
        h    = 20.0,
        B    = 1500.0,
        κ    = 2.0,
        ρ    = 3000.0,
        c    = 900.0
    )

    @time dry.solve(model;
        T     = 300.0,
        τ     = 2.0,
        t     = 2400.0,
        α     = 0.1,
        ε     = ε,
        M     = 20
    )

	residuals = SimulationResiduals(model.res[])

    fig1, ax1, p1 = plotsimulationresiduals(residuals; ε = ε)
	ax1.xticks = 0:3:12
	ax1.yticks = -11.2:0.2:-9.8
	xlims!(ax1, extrema(ax1.xticks.val))
	ylims!(ax1, extrema(ax1.yticks.val))
	
    fig2 = plotspheretemperature(model.grid.r, model.problem.x, model.B(1000000))
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
# ╟─795bb843-1e8a-4151-a289-532f6282ada3
# ╟─294e84f7-d230-4de6-825b-b5a32068c2dd
# ╟─0c1e97e7-5293-4504-950f-0635d62b6092
# ╟─86ad23ee-f11d-45cd-8ef3-30de6afa0568
# ╠═f5f89e4e-8dd8-4269-8405-4e8833cfab24
# ╟─20a61a1b-0e71-4b8f-9687-d7e836a4831d
# ╟─1552e5db-4a33-47ca-a66f-fe4fafb40945
# ╟─883a0bef-9743-4dfd-8e63-aec333e6a5d5
# ╟─09ac1638-83a4-4b86-96ba-12d2ddd00ba6
# ╟─b7781285-0b98-439f-85b4-d1b9cf72919a
# ╟─1d74deee-9948-443c-a04a-7b1f81f53bd1
