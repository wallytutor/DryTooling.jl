### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 52d557a4-8ef1-4b50-a145-7995745502c6
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using Symbolics
    using PlutoUI
    include("parameters.jl")
end

# ╔═╡ f175ae20-641b-11ee-2bf8-e72783c8a834
md"""
# Porosity calculations

In what follows we develop a simple approach for estimating quantities required for the computation of quantities relevant to heat and momentum transfer in the counter-current porous reactor model. Expressions are manipulated symbolically and then numerical results are obtained.

$(TableOfContents())
"""

# ╔═╡ a677c786-874c-4cc7-9871-6016064797a3
md"""
## Geometrical considerations

We start by declaring the relevant symbols:

|     |    |
|:---:|:---|
``ϕ`` | Volumetric porosity fraction
``l`` | Characteristic size of blocks
``D`` | Depth of reactor cross-section
``W`` | Width of reactor cross-section
``H`` | Height vertical reactor
"""

# ╔═╡ 5ef67a89-2d81-4863-b16f-a0f0a3c6321f
args = @variables ϕ l D W H

# ╔═╡ 13daafac-abfa-476d-b845-a78cc864559e
md"""
The surface area ``Aᵦ`` and volume ``Vᵦ`` of a single cubic block is computed as:
"""

# ╔═╡ e48a03d9-7189-4ff0-95cf-6bb942cab369
Aᵦ = 6l^2

# ╔═╡ 915f13a0-c078-47e8-ba4e-2eb0ab8d7de8
Vᵦ = l^3

# ╔═╡ 1eabfb91-5d33-4f4a-a9b4-6f4def62a1d6
md"""
The cross section of the reactor is ``Aᵥ``, total volume ``Vᵥ``, and perimeter ``Pᵥ`` are:
"""

# ╔═╡ 4dfcfd32-3b16-48b1-81e4-4371bb11a731
Aᵥ = D * W

# ╔═╡ bd7ff6f9-1578-4780-8666-6ca6cfa0ede8
Vᵥ = H * Aᵥ

# ╔═╡ 1d886b01-f1c5-4d85-ad61-970ed562c0e8
Pᵥ = 2 * (D + W)

# ╔═╡ a68a22cc-8a2f-4bee-82d7-14e1182381ca
md"""
From porosity ``ϕ`` the total volume of solids inside the reactor is ``Vₛ`` is:
"""

# ╔═╡ d0e39584-6c13-4804-8adc-5b5eb962940f
Vₛ = (1 - ϕ) * Vᵥ

# ╔═╡ 9876504d-9732-47c1-9ea5-003e6a916c8f
md"""
The number of blocks ``Nₛ`` is estimated as the ratio between ``Vₛ`` and a block volume ``Vᵦ``
"""

# ╔═╡ a87f0302-8e77-4dea-9e10-0469ecfdf6e3
Nₛ = Vₛ / Vᵦ

# ╔═╡ 19ae6a86-684d-4dca-b75f-27b9bfe71b9c
md"""
Thus the surface area of gas-blocks can be estimated as ``Aₛ``:
"""

# ╔═╡ d003c1fb-2c61-4031-b7bc-9c11d6e2dcc2
Aₛ = Nₛ * Aᵦ

# ╔═╡ 5ff1e62f-2a34-4dee-b2de-12fa670212c1
md"""
The blocks perimeter *per unit height* required in the PFR model is simply
"""

# ╔═╡ 48d89f5c-80e7-4c79-869a-c9dcd87b5c77
Pₛ = Aₛ / H

# ╔═╡ 855ddbc6-cbb5-4e6f-a9ac-f7f2bfb3a8a3
md"""
For the model implementation the ratio of perimeters ``σ`` can be interesting
"""

# ╔═╡ ba293489-93fb-4e65-8f5b-446c0f13552d
σ = Pₛ / Pᵥ

# ╔═╡ 4dcf432f-a13e-4c25-8f30-5ac16befaf10
md"""
We conclude this step by compiling a function for oerimeter evaluation.
"""

# ╔═╡ 078c50a8-4efa-463d-9110-61a6a38431bd
"Evaluates the perimeter of solids per unit height of reactor."
relativeperimeter = build_function(σ, args...; expression = Val{false})

# ╔═╡ 5fcd5364-5d57-4711-9c63-8be40607aa56
md"""
## Statistical considerations
"""

# ╔═╡ bf57201f-ff60-408c-a15b-b4388c07da48


# ╔═╡ 0d81315f-9962-450e-b42a-af7525754408


# ╔═╡ 9b99bd14-8b19-4069-bac9-ebe147b6b515


# ╔═╡ c7eaa50c-0aaf-49c0-bd71-69246cf42bff


# ╔═╡ b4bfb25d-0acc-4411-b461-366cac5ee7f7


# ╔═╡ 975eddf5-9ced-47a8-8dff-8f90593c9feb


# ╔═╡ 5cb8cc6e-1097-47a9-a0c0-ec93fbf711ef
let
    σₙ = relativeperimeter(ϕₛ, blocksize, REACTOR.D, REACTOR.W, REACTOR.H)
    σₙ, σₙ * REACTOR.P
end

# ╔═╡ 21ab6e10-9438-46b1-8ac7-3d93e4c70d89
md"""
## Tools
"""

# ╔═╡ Cell order:
# ╟─f175ae20-641b-11ee-2bf8-e72783c8a834
# ╟─a677c786-874c-4cc7-9871-6016064797a3
# ╟─5ef67a89-2d81-4863-b16f-a0f0a3c6321f
# ╟─13daafac-abfa-476d-b845-a78cc864559e
# ╟─e48a03d9-7189-4ff0-95cf-6bb942cab369
# ╟─915f13a0-c078-47e8-ba4e-2eb0ab8d7de8
# ╟─1eabfb91-5d33-4f4a-a9b4-6f4def62a1d6
# ╟─4dfcfd32-3b16-48b1-81e4-4371bb11a731
# ╟─bd7ff6f9-1578-4780-8666-6ca6cfa0ede8
# ╟─1d886b01-f1c5-4d85-ad61-970ed562c0e8
# ╟─a68a22cc-8a2f-4bee-82d7-14e1182381ca
# ╟─d0e39584-6c13-4804-8adc-5b5eb962940f
# ╟─9876504d-9732-47c1-9ea5-003e6a916c8f
# ╟─a87f0302-8e77-4dea-9e10-0469ecfdf6e3
# ╟─19ae6a86-684d-4dca-b75f-27b9bfe71b9c
# ╟─d003c1fb-2c61-4031-b7bc-9c11d6e2dcc2
# ╟─5ff1e62f-2a34-4dee-b2de-12fa670212c1
# ╟─48d89f5c-80e7-4c79-869a-c9dcd87b5c77
# ╟─855ddbc6-cbb5-4e6f-a9ac-f7f2bfb3a8a3
# ╟─ba293489-93fb-4e65-8f5b-446c0f13552d
# ╟─4dcf432f-a13e-4c25-8f30-5ac16befaf10
# ╟─078c50a8-4efa-463d-9110-61a6a38431bd
# ╟─5fcd5364-5d57-4711-9c63-8be40607aa56
# ╠═bf57201f-ff60-408c-a15b-b4388c07da48
# ╠═0d81315f-9962-450e-b42a-af7525754408
# ╠═9b99bd14-8b19-4069-bac9-ebe147b6b515
# ╠═c7eaa50c-0aaf-49c0-bd71-69246cf42bff
# ╠═b4bfb25d-0acc-4411-b461-366cac5ee7f7
# ╠═975eddf5-9ced-47a8-8dff-8f90593c9feb
# ╠═5cb8cc6e-1097-47a9-a0c0-ec93fbf711ef
# ╟─21ab6e10-9438-46b1-8ac7-3d93e4c70d89
# ╠═52d557a4-8ef1-4b50-a145-7995745502c6
