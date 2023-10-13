# -*- coding: utf-8 -*-

""" Base type for linear algebra problems. """
abstract type AbstractMatrixProblem end

""" Base type for physical models. """
abstract type AbstractPhysicalModel end

""" Base type for transport models. """
abstract type AbstractTransportModel end

""" Base type for thermodynamic models. """
abstract type AbstractGasThermo end

""" Base type for thermodynamic models. """
abstract type AbstractSolidThermo end

""" Base type for transport models. """
abstract type AbstractSolidTransport end

""" Base type for solid materials. """
abstract type AbstractSolidMaterial end

""" Base type for simplified mixture substances. """
abstract type AbstractMixtureSubstance end

""" Base type for simplified mixture phases. """
abstract type AbstractMixturePhase end
