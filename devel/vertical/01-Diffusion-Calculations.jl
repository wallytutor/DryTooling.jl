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
"""

# ╔═╡ 9b5b60ef-a62b-400e-8726-ade2d01b4fe6


# ╔═╡ 049cd654-b64d-4226-95b6-ae99b75e2fc9


# ╔═╡ b3e19d01-7886-4ccc-89bb-eb5531ffb391


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
# ╠═9b5b60ef-a62b-400e-8726-ade2d01b4fe6
# ╠═049cd654-b64d-4226-95b6-ae99b75e2fc9
# ╠═b3e19d01-7886-4ccc-89bb-eb5531ffb391
# ╠═ee1427c2-6d9e-4d94-8ded-711c4d928257
# ╠═63d7a037-0523-4e39-b443-7cf73dc401af
# ╠═6d2bc969-3dcf-433b-8263-d560eab70638
# ╟─20a61a1b-0e71-4b8f-9687-d7e836a4831d
# ╟─1552e5db-4a33-47ca-a66f-fe4fafb40945
# ╟─b7781285-0b98-439f-85b4-d1b9cf72919a
