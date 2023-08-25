var documenterSearchIndex = {"docs":
[{"location":"DocsCanteraAPI/#CCantera-Cantera-C-API-for-Julia","page":"Cantera API","title":"CCantera - Cantera C-API for Julia","text":"","category":"section"},{"location":"DocsCanteraAPI/","page":"Cantera API","title":"Cantera API","text":"This is an experimental interface to Cantera library based on its C-API under beta testing version 3.0. This interface is at its early days and has developped for Julia >= 1.9.0 under Windows 10/11. When it is stable enough it will be published as a package and tested under other platforms.","category":"page"},{"location":"DocsCanteraAPI/#Useful-links","page":"Cantera API","title":"Useful links","text":"","category":"section"},{"location":"DocsCanteraAPI/","page":"Cantera API","title":"Cantera API","text":"Source files\nHeader files\nct.h\nctfunc.h\nctmultiphase.h\nctonedim.h\nctreactor.h\nctrpath.h\nctsurf.h","category":"page"},{"location":"DocsCanteraAPI/#Implementation-status","page":"Cantera API","title":"Implementation status","text":"","category":"section"},{"location":"DocsCanteraAPI/","page":"Cantera API","title":"Cantera API","text":"Status Header Function Module\nTested ct.h ct_appdelete inlined\nStruct ct.h soln_newSolution wrapped\n ct.h soln_newInterface waitlist\nTested ct.h soln_del inlined\nTested ct.h soln_name inlined\nStruct ct.h soln_thermo inlined\nStruct ct.h soln_kinetics inlined\nStruct ct.h soln_transport inlined\nTested ct.h soln_setTransportModel wrapped\nTested ct.h soln_nAdjacent inlined\n ct.h soln_adjacent waitlist\nTested ct.h thermo_newFromFile wrapped\nTested ct.h thermo_del inlined\nStruct ct.h thermo_nElements inlined\nStruct ct.h thermo_nSpecies inlined\nTested ct.h thermo_temperature inlined\nStruct ct.h thermo_setTemperature inlined\nTested ct.h thermo_density inlined\nTested ct.h thermo_setDensity inlined\nTested ct.h thermo_molarDensity inlined\nTested ct.h thermo_setMolarDensity inlined\nTested ct.h thermo_meanMolecularWeight inlined\nTested ct.h thermo_moleFraction inlined\nTested ct.h thermo_massFraction inlined\nStruct ct.h thermo_getMoleFractions inlined\nTested ct.h thermo_getMassFractions inlined\nStruct ct.h thermo_setMoleFractions inlined\nTested ct.h thermo_setMassFractions inlined\n ct.h thermo_setMoleFractionsByName \n ct.h thermo_setMassFractionsByName \n ct.h thermo_getAtomicWeights \n ct.h thermo_getMolecularWeights \n ct.h thermo_getCharges \n ct.h thermo_getElementName \n ct.h thermo_getSpeciesName \n ct.h thermo_getName \n ct.h thermo_setName \n ct.h thermo_elementIndex \n ct.h thermo_speciesIndex \n ct.h thermo_report \nTested ct.h thermo_print \n ct.h thermo_nAtoms \n ct.h thermo_addElement \n ct.h thermo_getEosType \nTo test ct.h thermo_refPressure inlined\nTo test ct.h thermo_minTemp inlined\nTo test ct.h thermo_maxTemp inlined\nTo test ct.h thermoenthalpymole inlined\nTo test ct.h thermointEnergymole inlined\nTo test ct.h thermoentropymole inlined\nTo test ct.h thermogibbsmole inlined\nTo test ct.h thermocpmole inlined\nTo test ct.h thermocvmole inlined\nTo test ct.h thermo_pressure inlined\nStruct ct.h thermo_setPressure inlined\nTo test ct.h thermoenthalpymass inlined\nTo test ct.h thermointEnergymass inlined\nTo test ct.h thermoentropymass inlined\nTo test ct.h thermogibbsmass inlined\nTo test ct.h thermocpmass inlined\nTo test ct.h thermocvmass inlined\nTo test ct.h thermo_electricPotential inlined\nTo test ct.h thermo_thermalExpansionCoeff inlined\nTo test ct.h thermo_isothermalCompressibility inlined\n ct.h thermo_chemPotentials \n ct.h thermogetEnthalpiesRT \n ct.h thermogetEntropiesR \n ct.h thermogetCpR \n ct.h thermo_setElectricPotential \n ct.h thermosetTP \n ct.h thermosetTD \n ct.h thermosetRP \n ct.h thermosetDP \n ct.h thermosetHP \n ct.h thermosetUV \n ct.h thermosetSV \n ct.h thermosetSP \n ct.h thermosetST \n ct.h thermosetTV \n ct.h thermosetPV \n ct.h thermosetUP \n ct.h thermosetVH \n ct.h thermosetTH \n ct.h thermosetSH \nTested ct.h thermo_equilibrate \nTo test ct.h thermo_critTemperature inlined\nTo test ct.h thermo_critPressure inlined\nTo test ct.h thermo_critDensity inlined\nTo test ct.h thermo_vaporFraction inlined\n ct.h thermo_satTemperature \n ct.h thermo_satPressure \n ct.h thermosetStatePsat \n ct.h thermosetStateTsat \n ct.h kin_newFromFile \nTo test ct.h kin_del inlined\nTo test ct.h kin_nSpecies inlined\nTo test ct.h kin_nReactions inlined\nTo test ct.h kin_nPhases inlined\n ct.h kin_phaseIndex \nTo test ct.h kin_reactionPhaseIndex inlined\n ct.h kin_reactantStoichCoeff \n ct.h kin_productStoichCoeff \n ct.h kin_getReactionType \n ct.h kin_getFwdRatesOfProgress \n ct.h kin_getRevRatesOfProgress \n ct.h kin_getNetRatesOfProgress \n ct.h kin_getEquilibriumConstants \n ct.h kin_getFwdRateConstants \n ct.h kin_getRevRateConstants \n ct.h kin_getDelta \n ct.h kin_getCreationRates \n ct.h kin_getDestructionRates \n ct.h kin_getNetProductionRates \n ct.h kin_getSourceTerms \nTo test ct.h kin_multiplier inlined\n ct.h kin_getReactionString \n ct.h kin_setMultiplier \nTo test ct.h kin_isReversible inlined\n ct.h kin_getType \nTo test ct.h kin_start inlined\n ct.h kin_speciesIndex \nTo test ct.h kin_advanceCoverages inlined\nTo test ct.h kin_phase inlined\nTo test ct.h trans_newDefault inlined\n ct.h trans_new \nTo test ct.h trans_del inlined\nTo test ct.h trans_transportModel inlined\nTo test ct.h trans_viscosity inlined\nTo test ct.h trans_electricalConductivity inlined\n ct.h trans_thermalConductivity \n ct.h trans_getThermalDiffCoeffs \n ct.h trans_getMixDiffCoeffs \n ct.h trans_getBinDiffCoeffs \n ct.h trans_getMultiDiffCoeffs \n ct.h trans_setParameters \n ct.h trans_getMolarFluxes \n ct.h trans_getMassFluxes \n ct.h ct_getCanteraError \n ct.h ct_setLogWriter \n ct.h ct_setLogCallback \n ct.h ct_addCanteraDirectory \n ct.h ct_getDataDirectories \n ct.h ct_getCanteraVersion \n ct.h ct_getGitCommit \nTested ct.h ctsuppressthermo_warnings inlined\nTested ct.h ctuselegacyrateconstants inlined\nTested ct.h ct_clearStorage inlined\nTested ct.h ct_resetStorage inlined\n ctfunc.h func_new \n ctfunc.h funcnewbasic \n ctfunc.h funcnewadvanced \n ctfunc.h funcnewcompound \n ctfunc.h funcnewmodified \nTo test ctfunc.h func_del inlined\n ctfunc.h func_type \nTo test ctfunc.h func_value inlined\nTo test ctfunc.h func_derivative inlined\nTo test ctfunc.h func_duplicate inlined\n ctfunc.h func_write \nTo test ctfunc.h ct_clearFunc inlined\nTo test ctmultiphase.h mix_new inlined\nTo test ctmultiphase.h mix_del inlined\nTo test ctmultiphase.h ct_clearMix inlined\n ctmultiphase.h mix_addPhase \nTo test ctmultiphase.h mix_init inlined\nTo test ctmultiphase.h mix_updatePhases inlined\nTo test ctmultiphase.h mix_nElements inlined\n ctmultiphase.h mix_elementIndex \n ctmultiphase.h mix_speciesIndex \nTo test ctmultiphase.h mix_nSpecies inlined\nTo test ctmultiphase.h mix_setTemperature inlined\nTo test ctmultiphase.h mix_temperature inlined\nTo test ctmultiphase.h mix_minTemp inlined\nTo test ctmultiphase.h mix_maxTemp inlined\nTo test ctmultiphase.h mix_charge inlined\nTo test ctmultiphase.h mix_phaseCharge inlined\nTo test ctmultiphase.h mix_setPressure inlined\nTo test ctmultiphase.h mix_pressure inlined\nTo test ctmultiphase.h mix_nAtoms inlined\nTo test ctmultiphase.h mix_nPhases inlined\nTo test ctmultiphase.h mix_phaseMoles inlined\n ctmultiphase.h mix_setPhaseMoles \n ctmultiphase.h mix_setMoles \n ctmultiphase.h mix_setMolesByName \nTo test ctmultiphase.h mix_speciesMoles inlined\nTo test ctmultiphase.h mix_elementMoles inlined\n ctmultiphase.h mix_equilibrate \n ctmultiphase.h mix_getChemPotentials \nTo test ctmultiphase.h mix_enthalpy inlined\nTo test ctmultiphase.h mix_entropy inlined\nTo test ctmultiphase.h mix_gibbs inlined\nTo test ctmultiphase.h mix_cp inlined\nTo test ctmultiphase.h mix_volume inlined\nTo test ctmultiphase.h mix_speciesPhaseIndex inlined\nTo test ctmultiphase.h mix_moleFraction inlined\nTo test ctonedim.h ct_clearOneDim inlined\n ctonedim.h domain_new \nTo test ctonedim.h domain_del inlined\nTo test ctonedim.h domain_type inlined\n ctonedim.h domain_type3 \nTo test ctonedim.h domain_index inlined\nTo test ctonedim.h domain_nComponents inlined\nTo test ctonedim.h domain_nPoints inlined\n ctonedim.h domain_componentName \n ctonedim.h domain_componentIndex \n ctonedim.h domain_setBounds \nTo test ctonedim.h domain_lowerBound inlined\nTo test ctonedim.h domain_upperBound inlined\n ctonedim.h domain_setSteadyTolerances \n ctonedim.h domain_setTransientTolerances \nTo test ctonedim.h domain_rtol inlined\nTo test ctonedim.h domain_atol inlined\n ctonedim.h domain_setupGrid \n ctonedim.h domain_setID \nTo test ctonedim.h domain_grid inlined\nTo test ctonedim.h bdry_setMdot inlined\nTo test ctonedim.h bdry_setTemperature inlined\nTo test ctonedim.h bdry_setSpreadRate inlined\n ctonedim.h bdry_setMoleFractions \nTo test ctonedim.h bdry_temperature inlined\nTo test ctonedim.h bdry_spreadRate inlined\nTo test ctonedim.h bdry_massFraction inlined\nTo test ctonedim.h bdry_mdot inlined\nTo test ctonedim.h reactingsurf_setkineticsmgr inlined\nTo test ctonedim.h reactingsurf_enableCoverageEqs inlined\nTo test ctonedim.h inlet_new inlined\nTo test ctonedim.h outlet_new inlined\nTo test ctonedim.h outletres_new inlined\nTo test ctonedim.h symm_new inlined\nTo test ctonedim.h surf_new inlined\nTo test ctonedim.h reactingsurf_new inlined\nTo test ctonedim.h inlet_setSpreadRate inlined\n ctonedim.h stflow_new \nTo test ctonedim.h stflow_setTransport inlined\nTo test ctonedim.h stflow_enableSoret inlined\nTo test ctonedim.h stflow_setPressure inlined\nTo test ctonedim.h stflow_pressure inlined\n ctonedim.h stflow_setFixedTempProfile \nTo test ctonedim.h stflow_solveEnergyEqn inlined\n ctonedim.h sim1D_new \nTo test ctonedim.h sim1D_del inlined\n ctonedim.h sim1D_setValue \n ctonedim.h sim1D_setProfile \n ctonedim.h sim1D_setFlatProfile \n ctonedim.h sim1D_show \n ctonedim.h sim1D_showSolution \n ctonedim.h sim1D_setTimeStep \nTo test ctonedim.h sim1D_getInitialSoln inlined\n ctonedim.h sim1D_solve \nTo test ctonedim.h sim1D_refine inlined\n ctonedim.h sim1D_setRefineCriteria \n ctonedim.h sim1D_setGridMin \n ctonedim.h sim1D_save \n ctonedim.h sim1D_restore \nTo test ctonedim.h sim1D_writeStats inlined\n ctonedim.h sim1D_domainIndex \n ctonedim.h sim1D_value \n ctonedim.h sim1D_workValue \nTo test ctonedim.h sim1D_eval inlined\nTo test ctonedim.h sim1D_setMaxJacAge inlined\nTo test ctonedim.h sim1D_setFixedTemperature inlined\n ctreactor.h reactor_new \nTo test ctreactor.h reactor_del inlined\nTo test ctreactor.h reactor_setInitialVolume inlined\nTo test ctreactor.h reactor_setChemistry inlined\nTo test ctreactor.h reactor_setEnergy inlined\nTo test ctreactor.h reactor_setThermoMgr inlined\nTo test ctreactor.h reactor_setKineticsMgr inlined\nTo test ctreactor.h reactor_insert inlined\nTo test ctreactor.h reactor_mass inlined\nTo test ctreactor.h reactor_volume inlined\nTo test ctreactor.h reactor_density inlined\nTo test ctreactor.h reactor_temperature inlined\nTo test ctreactor.h reactorenthalpymass inlined\nTo test ctreactor.h reactorintEnergymass inlined\nTo test ctreactor.h reactor_pressure inlined\nTo test ctreactor.h reactor_massFraction inlined\nTo test ctreactor.h reactor_nSensParams inlined\nTo test ctreactor.h reactor_addSensitivityReaction inlined\nTo test ctreactor.h flowReactor_setMassFlowRate inlined\nTo test ctreactor.h reactornet_new inlined\nTo test ctreactor.h reactornet_del inlined\nTo test ctreactor.h reactornet_setInitialTime inlined\nTo test ctreactor.h reactornet_setMaxTimeStep inlined\nTo test ctreactor.h reactornet_setTolerances inlined\nTo test ctreactor.h reactornet_setSensitivityTolerances inlined\nTo test ctreactor.h reactornet_addreactor inlined\nTo test ctreactor.h reactornet_advance inlined\nTo test ctreactor.h reactornet_step inlined\nTo test ctreactor.h reactornet_time inlined\nTo test ctreactor.h reactornet_rtol inlined\nTo test ctreactor.h reactornet_atol inlined\n ctreactor.h reactornet_sensitivity \n ctreactor.h flowdev_new \nTo test ctreactor.h flowdev_del inlined\nTo test ctreactor.h flowdev_install inlined\nTo test ctreactor.h flowdev_setMaster inlined\nTo test ctreactor.h flowdev_setPrimary inlined\nTo test ctreactor.h flowdev_massFlowRate inlined\nTo test ctreactor.h flowdev_setMassFlowCoeff inlined\nTo test ctreactor.h flowdev_setValveCoeff inlined\nTo test ctreactor.h flowdev_setPressureCoeff inlined\nTo test ctreactor.h flowdev_setPressureFunction inlined\nTo test ctreactor.h flowdev_setTimeFunction inlined\n ctreactor.h wall_new \nTo test ctreactor.h wall_del inlined\nTo test ctreactor.h wall_install inlined\nTo test ctreactor.h wall_vdot inlined\nTo test ctreactor.h wall_expansionRate inlined\nTo test ctreactor.h wall_Q inlined\nTo test ctreactor.h wall_heatRate inlined\nTo test ctreactor.h wall_area inlined\nTo test ctreactor.h wall_setArea inlined\nTo test ctreactor.h wall_setThermalResistance inlined\nTo test ctreactor.h wall_setHeatTransferCoeff inlined\nTo test ctreactor.h wall_setHeatFlux inlined\nTo test ctreactor.h wall_setExpansionRateCoeff inlined\nTo test ctreactor.h wall_setVelocity inlined\nTo test ctreactor.h wall_setEmissivity inlined\nTo test ctreactor.h wall_ready inlined\nTo test ctreactor.h reactorsurface_new inlined\nTo test ctreactor.h reactorsurface_del inlined\nTo test ctreactor.h reactorsurface_install inlined\nTo test ctreactor.h reactorsurface_setkinetics inlined\nTo test ctreactor.h reactorsurface_area inlined\nTo test ctreactor.h reactorsurface_setArea inlined\nTo test ctreactor.h reactorsurface_addSensitivityReaction inlined\nTo test ctreactor.h ct_clearReactors inlined\nTo test ctrpath.h rdiag_new inlined\nTo test ctrpath.h rdiag_del inlined\nTo test ctrpath.h rdiag_detailed inlined\nTo test ctrpath.h rdiag_brief inlined\nTo test ctrpath.h rdiag_setThreshold inlined\n ctrpath.h rdiag_setBoldColor \n ctrpath.h rdiag_setNormalColor \n ctrpath.h rdiag_setDashedColor \n ctrpath.h rdiag_setDotOptions \nTo test ctrpath.h rdiag_setBoldThreshold inlined\nTo test ctrpath.h rdiag_setNormalThreshold inlined\nTo test ctrpath.h rdiag_setLabelThreshold inlined\nTo test ctrpath.h rdiag_setScale inlined\nTo test ctrpath.h rdiag_setFlowType inlined\nTo test ctrpath.h rdiag_setArrowWidth inlined\n ctrpath.h rdiag_setTitle \n ctrpath.h rdiag_write \nTo test ctrpath.h rdiag_add inlined\n ctrpath.h rdiag_findMajor \n ctrpath.h rdiag_setFont \nTo test ctrpath.h rdiag_displayOnly inlined\nTo test ctrpath.h rbuild_new inlined\nTo test ctrpath.h rbuild_del inlined\n ctrpath.h rbuild_init \n ctrpath.h rbuild_build \nTo test ctrpath.h ct_clearReactionPath inlined\n ctsurf.h surf_setCoverages \n ctsurf.h surf_getCoverages \n ctsurf.h surf_setConcentrations \n ctsurf.h surf_getConcentrations \nTo test ctsurf.h surf_setSiteDensity inlined\nTo test ctsurf.h surf_siteDensity inlined\n ctsurf.h surf_setCoveragesByName ","category":"page"},{"location":"#DryTooling","page":"Home","title":"DryTooling","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DryTooling.","category":"page"},{"location":"","page":"Home","title":"Home","text":"I am often faced with using the same approach for different engineering and scientific problems, but I don't like repeating the same task again and again. This is where DryTooling.jl comes in. By adopting some principles of DRY in Julia, to a larger extent than its definition, it packages together models and workflows that are not available or validated elsewhere - and in some cases adapts existing models. The tools will progressively cover a broad range of numerical applications and data treatment, this package is in its early days from the migration of my old Python scripts and packages.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Also dry tooling is my favorite sport!","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See CITATION.bib for the relevant reference(s).","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DryTooling","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [\n    DryTooling\n]","category":"page"},{"location":"#DryTooling.GAS_CONSTANT","page":"Home","title":"DryTooling.GAS_CONSTANT","text":"Ideal gas constant [J/(mol.K)]. \n\n\n\n\n\n","category":"constant"},{"location":"#DryTooling.ONE_ATM","page":"Home","title":"DryTooling.ONE_ATM","text":"Atmospheric pressure at sea level [Pa]. \n\n\n\n\n\n","category":"constant"},{"location":"#DryTooling.STABLE_ELEMENTS_TABLE","page":"Home","title":"DryTooling.STABLE_ELEMENTS_TABLE","text":"Instantiation of stable elements table. \n\n\n\n\n\n","category":"constant"},{"location":"#DryTooling.TRANSPORT_MODELS","page":"Home","title":"DryTooling.TRANSPORT_MODELS","text":"Named access to transport models. \n\n\n\n\n\n","category":"constant"},{"location":"#DryTooling.ZERO_CELSIUS","page":"Home","title":"DryTooling.ZERO_CELSIUS","text":"Zero degrees Celsius in Kelvin for conversion [K]. \n\n\n\n\n\n","category":"constant"},{"location":"#DryTooling.AbstractGasThermo","page":"Home","title":"DryTooling.AbstractGasThermo","text":"Base type for thermodynamic models. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.AbstractTransportModel","page":"Home","title":"DryTooling.AbstractTransportModel","text":"Base type for transport models. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.ElementData","page":"Home","title":"DryTooling.ElementData","text":"Represents a chemical element. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.IdealGasMixture","page":"Home","title":"DryTooling.IdealGasMixture","text":"Ideal gas phase mixture model. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.IdealGasSpecies","page":"Home","title":"DryTooling.IdealGasSpecies","text":"Ideal gas phase species model. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.IdealGasThermo","page":"Home","title":"DryTooling.IdealGasThermo","text":"Ideal gas phase thermodynamics model. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.LennardJonesTransport","page":"Home","title":"DryTooling.LennardJonesTransport","text":"Lennard-Jones ideal gas transport model. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.RotaryKilnBedSolution","page":"Home","title":"DryTooling.RotaryKilnBedSolution","text":"RotaryKilnBedSolution\n\nDescription of a rotary kiln bed geometry computed from the solution of bed height along the kiln length. The main goal of the quantities computed here is their use with heat and mass transfer models for the simulation of rotary kiln process.\n\nInternal elements are initialized through the following constructor:\n\nRotaryKilnBedSolution(\n    z::Vector{Float64},\n    h::Vector{Float64},\n    R::Float64,\n    Φ::Float64\n)\n\nWhere parameters are given as:\n\n- `z`: solution coordinates over length, [m].\n- `h`: bed profile solution over length, [m].\n- `R`: kiln internal radius, [m].\n- `Φ`: kiln feed rate, [m³/s].\n\nz::Vector{Float64}: Solution coordinates [m]\nh::Vector{Float64}: Solution bed height [m]\nθ::Vector{Float64}: View angle from kiln center [rad]\nl::Vector{Float64}: Bed-freeboard cord length [m]\nA::Vector{Float64}: Local bed cross section area [m²]\nη::Vector{Float64}: Local loading based on height [-]\nηₘ::Float64: Mean loading of kiln [%]\nV::Float64: Bed integral volume [m³]\nτ::Float64: Residence time of particles\n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.StableElementsTable","page":"Home","title":"DryTooling.StableElementsTable","text":"Periodic table of stable chemical elements. \n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.SymbolicLinearKramersModel","page":"Home","title":"DryTooling.SymbolicLinearKramersModel","text":"SymbolicLinearKramersModel\n\nCreates a reusable linear Kramers model for rotary kiln simulation.\n\nImplements the ordinary differential equation for prediction of bed height profile in a rotary kiln as proposed by Kramers and Croockewite (1952) [1]. Its goal is to be used as a process support tool or to integrate more complex models requiring integration of the bed profile.\n\n[1]: Kramers et al., 1952\n\nR::Symbolics.Num: Symbolic kiln internal radius\nΦ::Symbolics.Num: Symbolic kiln feed rate\nω::Symbolics.Num: Symbolic kiln rotation rate\nβ::Symbolics.Num: Symbolic kiln slope\nγ::Symbolics.Num: Symbolic solids repose angle\nz::Symbolics.Num: Symbolic kiln axial coordinates\nh::Symbolics.Num: Symbolic bed height profile\nsys::ModelingToolkit.ODESystem: Problem ordinary differential equation\n\n\n\n\n\n","category":"type"},{"location":"#DryTooling.densitymass-Tuple{DryTooling.IdealGasSolution}","page":"Home","title":"DryTooling.densitymass","text":"Mixture specific mass [kg/m³]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.element-Tuple{String}","page":"Home","title":"DryTooling.element","text":"Retrieve an element by name. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.elementmass-Tuple{String}","page":"Home","title":"DryTooling.elementmass","text":"Retrieve atomic mass of element from atomic symbol [g/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.enthalpymass-Tuple{DryTooling.IdealGasSpecies, Any}","page":"Home","title":"DryTooling.enthalpymass","text":"Species enthalpy in mass units [J/kg]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.enthalpymole-Tuple{DryTooling.IdealGasSpecies, Any}","page":"Home","title":"DryTooling.enthalpymole","text":"Species enthalpy in mole units [J/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.getnameditem-Tuple{Any, Any}","page":"Home","title":"DryTooling.getnameditem","text":"Query first item matching name in dictionary. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.getthermo-NTuple{6, Any}","page":"Home","title":"DryTooling.getthermo","text":"Create specific heat and enthalpy functions for species. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.heaviside-Tuple{Any}","page":"Home","title":"DryTooling.heaviside","text":"Automatic differentiable Heaviside function. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.interval-Tuple{Any}","page":"Home","title":"DryTooling.interval","text":"Returns 1 if x  (a b), 1/2 for x = a  x = b, or 0 . \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.makestepwise1d-Tuple{Any, Any, Any}","page":"Home","title":"DryTooling.makestepwise1d","text":"makestepwise1d(lo, hi, xc)\n\nCreates an univariate function that is composed of two parts, the first evaluated before a critical domain point xc, and the seconda above that value. This is often required, for instance, for the evaluation of NASA polynomials for thermodynamic properties. If differentiable, then the returned function is compatible with symbolic argument as required when using package ModelingToolkit, etc.\n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.mass-Tuple{DryTooling.ElementData}","page":"Home","title":"DryTooling.mass","text":"Retrieve atomic mass of element [g/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.mass-Tuple{DryTooling.IdealGasSpecies}","page":"Home","title":"DryTooling.mass","text":"Retrieve atomic mass of species [kg/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.massfraction2molefraction-Tuple{Vector{Float64}, Any}","page":"Home","title":"DryTooling.massfraction2molefraction","text":"Convert mass fractions to mole fractions. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.massfractions-Tuple{DryTooling.IdealGasSolution}","page":"Home","title":"DryTooling.massfractions","text":"Mixture composition in mole fractions. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.meanmolecularmass-Tuple{DryTooling.IdealGasSolution}","page":"Home","title":"DryTooling.meanmolecularmass","text":"Mixture mean molecular mass [kg/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.meanmolecularmass-Tuple{Vector{Float64}, Any}","page":"Home","title":"DryTooling.meanmolecularmass","text":"Mixture mean molecular mass [kg/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.molefraction2massfraction-Tuple{Vector{Float64}, Any}","page":"Home","title":"DryTooling.molefraction2massfraction","text":"Convert mole fractions to mass fractions. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.molefractions-Tuple{DryTooling.IdealGasSolution}","page":"Home","title":"DryTooling.molefractions","text":"Mixture composition in mole fractions. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.nasa7enthapy-Tuple{Any, Any}","page":"Home","title":"DryTooling.nasa7enthapy","text":"Molar enthalpy from NASA7 polynomial [J/mol]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.nasa7specificheat-Tuple{Any, Any}","page":"Home","title":"DryTooling.nasa7specificheat","text":"Molar specific heat from NASA7 polynomial [J/(mol.K)]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.plotlinearkramersmodel-Tuple{RotaryKilnBedSolution}","page":"Home","title":"DryTooling.plotlinearkramersmodel","text":"plotlinearkramersmodel(\n    model::SolutionLinearKramersModel;\n    normz::Bool = false,\n    normh::Bool = false\n)::Any\n\nDisplay plot of model solution for rotary kiln bed profile. Arguments normz and normh control whether z-coordinate and bed height must be normalized, respectively.\n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.solvelinearkramersmodel-Tuple{}","page":"Home","title":"DryTooling.solvelinearkramersmodel","text":"solvelinearkramersmodel(;\n    model::SymbolicLinearKramersModel,\n    L::Float64,\n    R::Float64,\n    Φ::Float64,\n    ω::Float64,\n    β::Float64,\n    γ::Float64,\n    d::Float64,\n    solver::Any = Tsit5(),\n    rtol::Float64 = 1.0e-08,\n    atol::Float64 = 1.0e-08\n)\n\nIntegrates an instance of SymbolicLinearKramersModel.\n\nImportant: inputs must be provided in international system (SI) units as a better physical practice. The only exception is the rotation rate ω provided in revolution multiples. If the discharge end is held by a dam, its height must be provided instead of the particle size, as it is used as the ODE initial condition.\n\nmodel: a symbolic kiln model.\nL: kiln length, [m].\nR: kiln internal radius, [m].\nΦ: kiln feed rate, [m³/s].\nω: kiln rotation rate, [rev/s].\nβ: kiln slope, [rad].\nγ: solids repose angle, [rad].\nd: particle size or dam height, [m].\n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.specificheatmass-Tuple{DryTooling.IdealGasSolution}","page":"Home","title":"DryTooling.specificheatmass","text":"Mixture mass-averaged specific heat [J/(kg.K)]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.specificheatmass-Tuple{DryTooling.IdealGasSpecies, Any}","page":"Home","title":"DryTooling.specificheatmass","text":"Species specific heat in mass units [J/(kg.K)]. \n\n\n\n\n\n","category":"method"},{"location":"#DryTooling.specificheatmole-Tuple{DryTooling.IdealGasSpecies, Any}","page":"Home","title":"DryTooling.specificheatmole","text":"Species specific heat in mole units [J/(mol.K)]. \n\n\n\n\n\n","category":"method"}]
}
