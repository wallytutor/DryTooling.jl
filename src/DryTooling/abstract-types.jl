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

""" Base type for """
abstract type AbstractSolidMaterial end

""" """
abstract type Substance end

""" """
abstract type Mixture end