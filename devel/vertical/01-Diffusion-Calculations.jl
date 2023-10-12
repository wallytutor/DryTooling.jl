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

    using Symbolics
    using PlutoUI
    include("parameters.jl")
end

# ╔═╡ b307aa70-641c-11ee-27cb-b151d31fffc9
md"""
# Diffusion Calculations

$(TableOfContents())
"""

# ╔═╡ bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
md"""
## Sphere heating model

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

To get a solution with finite volume method we start with an integration over the derivative variables. In what follows an implicit time-stepping method is assumed without further discussions.

```math
\int_{r_s}^{r_n}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=
\int_{r_s}^{r_n}\int_{0}^{\tau}
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)
```
"""

# ╔═╡ 9b5b60ef-a62b-400e-8726-ade2d01b4fe6
md"""
```math
\int_{V}\int_{t_0}^{t_1}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}dtdV=
\int_{V}\int_{t_0}^{t_1}
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)dtdV
```

```math
dV=r^2\sin\phi{}dr{}d\theta{}d\phi
```

```math
\int_{0}^{\pi}\int_{0}^{2\pi}\sin\phi{}d\theta{}d\phi=4\pi
```

```math
\int_{r_s}^{r_n}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}r^2dtdr=
\int_{r_s}^{r_n}\int_{0}^{\tau}
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)r^2dtdr
```
"""

# ╔═╡ 049cd654-b64d-4226-95b6-ae99b75e2fc9
md"""
```math
\rho{}c_{p}\left[T_P^{\tau}-T_P^{0}\right]
\left(\frac{r^3}{3}\right)\biggr\vert_{r_s}^{r_n}=
\tau
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)\biggr\vert_{r_s}^{r_n}
```


```math
\begin{align}
\alpha_{sn} & = \frac{\rho{}c_{p}}{3\tau}\left(r_n^3-r_s^3\right)\\[8pt]
\beta_{k}   & = r_k^2k_k\\[8pt]
\delta_{K}  & = \vert{}r_{K}-r_{P}\vert
\end{align}
```

```math
\alpha_{sn}T_P^{\tau}-\alpha_{sn}T_P^{0}=
\frac{\beta_{n}}{\delta_{N}}\left(T_N-T_P\right)-
\frac{\beta_{s}}{\delta_{S}}\left(T_P-T_S\right)
```
"""

# ╔═╡ b3e19d01-7886-4ccc-89bb-eb5531ffb391
md"""
```math
-\frac{\beta_{s}}{\delta_{S}}T_S
+\left(\alpha_{sn}+\frac{\beta_{n}}{\delta_{N}}+\frac{\beta_{s}}{\delta_{S}}\right)T_P-
\frac{\beta_{n}}{\delta_{N}}T_N=
\alpha_{sn}T_P^{0}
```

```math
\begin{align}
a_{S} & = -\frac{\beta_{s}}{\delta_{S}}\\[8pt]
a_{N} & = -\frac{\beta_{n}}{\delta_{N}}\\[8pt]
a_{P} & = \alpha_{sn}+\frac{\beta_{n}}{\delta_{N}}+\frac{\beta_{s}}{\delta_{S}}
\end{align}
```

```math
a_{S}T_S + a_{P}T_P + a_{N}T_N = \alpha_{sn}T_P^{0}
```
"""

# ╔═╡ ee1427c2-6d9e-4d94-8ded-711c4d928257


# ╔═╡ 63d7a037-0523-4e39-b443-7cf73dc401af


# ╔═╡ 6d2bc969-3dcf-433b-8263-d560eab70638


# ╔═╡ 20a61a1b-0e71-4b8f-9687-d7e836a4831d
md"""
## Advection-diffusion plug-flow
"""

# ╔═╡ 1552e5db-4a33-47ca-a66f-fe4fafb40945
md"""
## Tools
"""

# ╔═╡ Cell order:
# ╟─b307aa70-641c-11ee-27cb-b151d31fffc9
# ╟─bd9279ee-3e02-4700-a3c8-a1fed7dd1d78
# ╟─9b5b60ef-a62b-400e-8726-ade2d01b4fe6
# ╟─049cd654-b64d-4226-95b6-ae99b75e2fc9
# ╟─b3e19d01-7886-4ccc-89bb-eb5531ffb391
# ╠═ee1427c2-6d9e-4d94-8ded-711c4d928257
# ╠═63d7a037-0523-4e39-b443-7cf73dc401af
# ╠═6d2bc969-3dcf-433b-8263-d560eab70638
# ╟─20a61a1b-0e71-4b8f-9687-d7e836a4831d
# ╟─1552e5db-4a33-47ca-a66f-fe4fafb40945
# ╟─b7781285-0b98-439f-85b4-d1b9cf72919a
