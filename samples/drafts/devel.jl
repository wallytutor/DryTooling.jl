### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 696971cb-3122-4bae-82b1-bbc87a413875
begin
	cantera = "C:\\Program Files\\Cantera\\bin\\cantera_shared.dll"
	
	ENV["CANTERA_SHARED"] = cantera
	@assert haskey(ENV, "CANTERA_SHARED")
	
	include("..\\src\\CCantera.jl")
end;

# ╔═╡ faae1361-88ce-4fe8-a012-ce3fc68d5264
md"""
# Interface com C

Neste tutorial vamos ilustrar como usar uma biblioteca de link dinâmico (DLL) em Julia. Essa atividade é recorrente quando desejamos realizar interface de código externo -- normalmente bibliotecas mais antigas ou aceitas como de alta performance -- com código desenvolvido para aplicações específicas.

Um caso que encontramos com frequência na temática global do livro texto é a avaliação de propriedades de transporte em produtos de combustão em termos de uma composição local. Isso é possível de ser realizado através do pacote Cantera, que pelo momento não é disponível em Julia [^1]. As interfaces em C da *DLL* [^2] de Cantera incluem o arquivo de cabeçalho [clib_defs.h](https://cantera.org/documentation/docs-2.6/doxygen/html/dd/d7b/clib__defs_8h.html).

[^1]: Funções de Cantera poderiam ser facilmente acessadas através da interface em Python, mas isso adicionaria a necessidade de se utilizar Conda e a interface de chamada de funções em Python. Por razões evidentes e para manter o sentido deste tutorial vamos negligenciar essa possibilidade.
[^2]: Nos sistemas operacionais Linux (e similares) e Mac OS a terminologia para uma biblioteca compartilhada utiliza a extensão `.so`. Para nosso deleito Julia leva em conta essas diferenças e veremos que o nome da biblioteca não necessita especificar a extensão do arquivo.
"""

# ╔═╡ d195cfb6-ea77-4103-b187-4d97844dd979
import .CCantera as ct;

# ╔═╡ f2cafce4-e13d-47f3-b903-b3ce3d881557
begin
	sol = ct.Solution("gri30.yaml", "gri30", "mixture-averaged")

	Xᵣ = zeros(sol.nspecies)
	Xᵣ[1] = 1.0

	Tᵣ = 3500.0
	Pᵣ = 50000.0
	
	ct.set_TPX!(sol, Tᵣ, Pᵣ, Xᵣ; norm = true)

	@assert ct.gettemperature(sol) ≈ Tᵣ
	@assert ct.getpressure(sol) ≈ Pᵣ
	@assert all(ct.getmolefractions(sol) ≈ Xᵣ)
end

# ╔═╡ 5d60edce-811c-462b-8e75-bfc626568b05
ct.getpressure(sol)

# ╔═╡ 9bf2ce93-9c8c-4ce7-bfca-43468ee086ea
# ct.soln_test();

# ╔═╡ 46a7c1c8-2bb6-4f10-b06f-0a8fc5d64726
# ct.thermo_test();

# ╔═╡ abf2d0ac-273f-4d79-8a86-15d6477fa8a7
# function soln_report(n; newline = false)
#     if newline
#         println("")
#     end

#     println(string(
#         "Solution $(n):",
#         "\n- name ........... $(soln_name(n))",
#         "\n- thermo ......... $(soln_thermo(n))",
#         "\n- kinetics ....... $(soln_kinetics(n))",
#         "\n- transport ...... $(soln_transport(n))",
#         "\n- nAdjacent ...... $(soln_nAdjacent(n))"
#     ))
# end

# function soln_test()
#     ct_appdelete();
#     sol = Solution("gri30.yaml", "gri30", "mixture-averaged")
#     n = sol.index.solution

#     soln_report(n)
#     soln_setTransportModel(n, "ionized-gas")
#     soln_report(n, newline = true)
#     soln_del(n)
# end

# function thermo_report(n; newline = false, fractions = false)
#     if newline
#         println("")
#     end

#     println(string(
#         "Thermo $(n):",
#         "\n- nSpecies ....... $(thermo_nSpecies(n))",
#         "\n- nElements ...... $(thermo_nElements(n))",
#         "\n- Temperature .... $(thermo_temperature(n))",
#         "\n- Density ........ $(thermo_density(n))",
#         "\n- Molar density .. $(thermo_molarDensity(n))",
#         "\n- Mol. weight .... $(thermo_meanMolecularWeight(n))",
#     ));

#     if fractions
#         println("");
#         nspecies = thermo_nSpecies(n);

#         # Set all fractions randomly and normalize.
#         mrnd = rand(nspecies);
#         thermo_setMoleFractions(n, nspecies, mrnd, 1);
#         thermo_setMassFractions(n, nspecies, mrnd, 1);

#         xarr = zeros(nspecies);
#         yarr = zeros(nspecies);

#         thermo_getMoleFractions(n, nspecies, xarr);
#         thermo_getMassFractions(n, nspecies, yarr);

#         for k in 0:nspecies-1;
#             xnum = thermo_moleFraction(n, k);
#             ynum = thermo_massFraction(n, k);
            
#             @assert xarr[k+1] ≈ xnum
#             @assert yarr[k+1] ≈ ynum

#             x = @sprintf("%.6e", xnum);
#             y = @sprintf("%.6e", ynum);
#             println(" ... $(x) | $(y)")
#         end
#     end
# end

# function thermo_test()
#     ct_appdelete();
#     n = thermo_newFromFile("gri30.yaml", "gri30");
#     thermo_report(n);

#     thermo_setTemperature(n, 400.0);
#     thermo_report(n, newline = true, fractions = false);

#     thermo_setDensity(n, 10.0);
#     thermo_report(n, newline = true, fractions = false);

#     thermo_setMolarDensity(n, 10.0);
#     thermo_report(n, newline = true, fractions = true);

#     thermo_del(n);
# end

# ╔═╡ 9f0e5fb9-c44f-4246-bf0d-de1ed58d2201
# Cantera.ct_addCanteraDirectory("C:")

# ╔═╡ 7de960af-bddf-4915-b1ac-b571e68b4b76
# Cantera.ct_getDataDirectories()

# ╔═╡ 5c3aa251-e68c-4444-a63a-42a756a01eca
# Cantera.ct_getCanteraVersion()

# ╔═╡ 7abea6f1-5e61-4e18-9ce8-834e4d5bb278
# Cantera.ct_getGitCommit()

# ╔═╡ a2a0fb1d-8db2-4695-ba11-d7773427cbc7
# Cantera.ct_suppress_thermo_warnings(1)

# ╔═╡ 2f15da09-12fc-4437-9a1b-06c7247b302d
# Cantera.ct_use_legacy_rate_constants(1)

# ╔═╡ 4365627a-1eed-493f-8d24-ab8e7523a3ae
# Cantera.ct_clearStorage()

# ╔═╡ 710ec0c4-2aaf-4525-b88c-457041296403
# function JuliaCanteraExample()

#     xml_file = ccall((:JL_xml_get_XML_File, CCantera),
#                      Int32, (Cstring, Int32),
#                      "gri30.xml", 0)
#     testReturn("JL_xml_get_XML_File", xml_file)


#     phase_node = ccall((:JL_xml_findID, CCantera),
#                       Int32, (Int32, Cstring),
#                       xml_file, "gri30_mix")
#     testReturn("JL_xml_findID", phase_node)


#     thermo = ccall((:JL_thermo_newFromXML, CCantera),
#                    Int32, (Int32,),
#                    phase_node)
#     testReturn("JL_thermo_newFromXML", thermo)


#     nsp = ccall((:JL_thermo_nSpecies, CCantera),
#                 Int32, (Int32,),
#                 thermo)
#     testReturn("JL_thermo_nSpecies", nsp, eq=true, value=53)


#     ret = ccall((:JL_thermo_setTemperature, CCantera),
#                 Int32, (Int32, Float64),
#                 thermo, 500.0)
#     testReturn("JL_thermo_setTemperature", ret, eq=true, value=0)


#     ret = ccall((:JL_thermo_setPressure, CCantera),
#                 Int32, (Int32, Float64),
#                 thermo, 5.0 * 101325.0)
#     testReturn("JL_thermo_setPressure", ret, eq=true, value=0)


#     ret = ccall((:JL_thermo_setMoleFractionsByName, CCantera),
#                 Int32, (Int32, Cstring),
#                 thermo, "CH4:1.0, O2:2.0, N2:7.52")
#     testReturn("JL_thermo_setMoleFractionsByName", ret, eq=true, value=0)


#     ret = ccall((:JL_thermo_equilibrate, CCantera),
#                 Int32, (Int32, Cstring, Int32, Float64, Int32, Int32, Int32),
#                 thermo, "HP", 0, 1e-9, 50000, 1000, 0);
#     testReturn("JL_thermo_equilibrate", ret, eq=true, value=0)


#     T = ccall((:JL_thermo_temperature, CCantera),
#               Float64, (Int32,),
#               thermo)
#     ret = T > 2200 && T < 2300
#     testReturn("JL_thermo_equilibrate", ret, eq=true, value=true)


#     ret = ccall((:JL_thermo_print, CCantera),
#               Int32, (Int32, Int32, Int32),
#               thermo, 1, 0)
#     testReturn("JL_thermo_print", ret, eq=true, value=0)


#     kin = ccall((:JL_kin_newFromXML, CCantera),
#               Int32, (Int32, Int32, Int32, Int32, Int32, Int32),
#               phase_node, thermo, 0, 0, 0, 0)
#     testReturn("JL_kin_newFromXML", thermo)


#     nr = ccall((:JL_kin_nReactions, CCantera),
#               UInt32, (Int32,),
#               kin)
#     testReturn("JL_kin_nReactions", nr, eq=true, value=325)


#     ret = ccall((:JL_thermo_setTemperature, CCantera),
#                 Int32, (Int32, Float64),
#                 thermo, T - 200.0)
#     testReturn("JL_thermo_setTemperature", ret, eq=true, value=0)

#     buf::String = " " ^ 50
#     ropf::Vector{Float64} = Vector{Float64}(325)

#     println("\n                   Reaction           Forward ROP\n")

#     # ATTENTION: need to declare Ptr{Void} for double*!
#     ret = ccall((:JL_kin_getFwdRatesOfProgress, CCantera),
#                 Int32, (Int32, UInt32, Ptr{Void}),
#                 kin, 325, ropf)

#     for n=1:325
#         # Here had to use Ptr{UInt8} instead of Cstring
#         ret = ccall((:JL_kin_getReactionString, CCantera),
#                     Int32, (Int32, Int32, Int32, Ptr{UInt8}),
#                     kin, n, sizeof(buf), buf)
#         rate = ropf[n]
#         #TODO Consider using Formatting
#         println("$buf $rate\n")
#     end


#     println("\n  Species    Mix diff coeff\n")
#     tran = ccall((:JL_trans_new, CCantera),
#                  Int32, (Cstring, Int32, Int32),
#                  "Mix", thermo, 0)

#     dkm::Vector{Float64} = Vector{Float64}(53)
#     ret = ccall((:JL_trans_getMixDiffCoeffs, CCantera),
#                  Int32, (Int32, Int32, Ptr{Void}),
#                  tran, 53, dkm)

#     #FIXME this should not be necessary!
#     buf = " " ^ 50

#     for k=1:nsp

#         ret = ccall((:JL_thermo_getSpeciesName, CCantera),
#                      Int32, (Int32, UInt32, UInt32, Ptr{UInt8}),
#                      thermo, k, sizeof(buf), buf)
#          frac = dkm[k]
#          #TODO Consider using Formatting
#          println("$buf $frac\n")
#     end


#     ret = ccall((:JL_thermo_setTemperature, CCantera),
#                 Int32, (Int32, Float64),
#                 thermo, 1050.0)
#     testReturn("JL_thermo_setTemperature", ret, eq=true, value=0)


#     ret = ccall((:JL_thermo_setPressure, CCantera),
#                 Int32, (Int32, Float64),
#                 thermo, 5.0 * 101325.0)
#     testReturn("JL_thermo_setPressure", ret, eq=true, value=0)


#     ret = ccall((:JL_thermo_setMoleFractionsByName, CCantera),
#                 Int32, (Int32, Cstring),
#                 thermo, "CH4:1.0, O2:2.0, N2:7.52")
#     testReturn("JL_thermo_setMoleFractionsByName", ret, eq=true, value=0)



#     reactor = ccall((:JL_reactor_new, CCantera),
#                     Int32, (Int32,),
#                     5)


#     net = ccall((:JL_reactornet_new, CCantera),
#                     Int32, ())


#     ret = ccall((:JL_reactor_setThermoMgr, CCantera),
#                 Int32, (Int32, Int32),
#                 reactor, thermo)
#     testReturn("JL_reactor_setThermoMgr", ret, eq=true, value=0)


#     ret = ccall((:JL_reactor_setKineticsMgr, CCantera),
#                 Int32, (Int32, Int32),
#                 reactor, kin)
#     testReturn("JL_reactor_setKineticsMgr", ret, eq=true, value=0)


#     ret = ccall((:JL_reactornet_addreactor, CCantera),
#                 Int32, (Int32, Int32),
#                 net, reactor)
#     testReturn("JL_reactornet_addreactor", ret, eq=true, value=0)


#     println("\ntime       Temperature\n");

#     t::Float64 = 0.0

#     while t < 0.1 && ret == 0
#         T = ccall((:JL_reactor_temperature, CCantera),
#                   Float64, (Int32,),
#                   reactor)

#         t = ccall((:JL_reactornet_time, CCantera),
#                   Float64, (Int32,),
#                   net)

#         println("$t  $T")


#         ret = ccall((:JL_reactornet_advance, CCantera),
#                     Int32, (Int32, Float64),
#                     net, t + 5.0e-03)
#         testReturn("JL_reactornet_advance", ret, eq=true, value=0)
#     end


#     ret = ccall((:JL_ct_appdelete, CCantera), Int32, ())
#     testReturn("JL_ct_appdelete", ret, eq=true, value=0)

# end

# ╔═╡ Cell order:
# ╟─faae1361-88ce-4fe8-a012-ce3fc68d5264
# ╠═696971cb-3122-4bae-82b1-bbc87a413875
# ╠═d195cfb6-ea77-4103-b187-4d97844dd979
# ╠═f2cafce4-e13d-47f3-b903-b3ce3d881557
# ╠═5d60edce-811c-462b-8e75-bfc626568b05
# ╠═9bf2ce93-9c8c-4ce7-bfca-43468ee086ea
# ╠═46a7c1c8-2bb6-4f10-b06f-0a8fc5d64726
# ╠═abf2d0ac-273f-4d79-8a86-15d6477fa8a7
# ╠═9f0e5fb9-c44f-4246-bf0d-de1ed58d2201
# ╠═7de960af-bddf-4915-b1ac-b571e68b4b76
# ╠═5c3aa251-e68c-4444-a63a-42a756a01eca
# ╠═7abea6f1-5e61-4e18-9ce8-834e4d5bb278
# ╠═a2a0fb1d-8db2-4695-ba11-d7773427cbc7
# ╠═2f15da09-12fc-4437-9a1b-06c7247b302d
# ╠═4365627a-1eed-493f-8d24-ab8e7523a3ae
# ╠═710ec0c4-2aaf-4525-b88c-457041296403
