# -*- coding: utf-8 -*-
# # Mixture properties
#
# This notebook aims at fitting simple polynomials to the thermophysical properties of combustion flue gases and their dilution in carbon dioxide or water issued from calcination processes. A final model of a 3-component mixture is then proposed and validated.

import json
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Both manometric pressure and temperature range will be shared across the mixtures.

P = 12500.0
T = np.linspace(400, 2800, 50)


# Implementation of *pure substance* properties fitting is provided by `SubstanceFit`.

class SubstanceFit:
    """ Fit properties of mixture in temperature range.

    Parameters
    ----------
    T: list[float]
        Temperature array for fitting properties [K].
    P: float
        Pressure with reference to `ct.one_atm` [Pa].
    Y: dict[str, float]
        Dictionary of species mass fractions [-].
    degs: dict[str, int]
        Order of polynomials to fit. Supported entries for properties
        are `cp_mass`, `viscosity` and `thermal_conductivity`.
    mech: str = "gri30.yaml"
        Kinetics mechanism for creating solution.
    """
    def __init__(
            self,
            T: list[float],
            P: float,
            Y: dict[str, float],
            degs: dict[str, int] = {
                "cp_mass": 4,
                "viscosity": 4,
                "thermal_conductivity": 4
            },
            mech: str = "gri30.yaml"
        ):
        gas = ct.Solution(mech)
        sol = ct.SolutionArray(gas, shape=T.shape)
        sol.TPY = T, ct.one_atm + P, Y

        self._sol = sol
        self._mw = sol.mean_molecular_weight[1]
        self._Y = Y

        def fitrefs(p, d):
            return np.polyfit(sol.T, self._refs[p], deg=d)

        self._refs = {p: getattr(sol, p) for p, d in degs.items()}
        self._fits = {p: fitrefs(p, d)   for p, d in degs.items()}
        self._coef = {
            "cp": list(reversed(self._fits["cp_mass"])),
            "kg": list(reversed(self._fits["thermal_conductivity"])),
            "mu": list(reversed(self._fits["viscosity"]))
        }

    def plotsubstancefit(
            self,
            dT: float = 400
        ):
        """ Graphic evaluation of fitting for a substance. """
        T = self._sol.T
        xlims = T[0], T[-1]
        xtick = np.arange(T[0], T[-1]+1, 400.0)

        cpr = self._refs["cp_mass"]
        mur = self._refs["viscosity"]
        kgr = self._refs["thermal_conductivity"]

        cpf = np.polyval(self._fits["cp_mass"], T)
        muf = np.polyval(self._fits["viscosity"], T)
        kgf = np.polyval(self._fits["thermal_conductivity"], T)

        plt.close("all")
        plt.style.use("seaborn-white")
        fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

        ax[0].plot(T, cpr, "k")
        ax[1].plot(T, mur, "k")
        ax[2].plot(T, kgr, "k")

        ax[0].plot(T, cpf, "r")
        ax[1].plot(T, muf, "r")
        ax[2].plot(T, kgf, "r")

        ax[0].grid(linestyle=":")
        ax[1].grid(linestyle=":")
        ax[2].grid(linestyle=":")

        ax[0].set_ylabel("Specific heat [J/(kg.K)]")
        ax[1].set_ylabel("Viscosity [Pa.s]")
        ax[2].set_ylabel("Conductivity [W/(m.K)]")

        ax[0].set_xlim(xlims)
        ax[0].set_xticks(xtick)

        fig.tight_layout()

    def to_solution(self, T=1000, P=P):
        """ Get simple solution from mixture. """
        sol = ct.Solution(self._sol.source)
        sol.TPY = T, P, self._Y
        return sol

    @property
    def coefficients(self):
        """ Polynomial coefficients in decreasing order. """
        return self._coef


# ## Carbon dioxide

mix0 = SubstanceFit(T, P, {"CO2": 1.0})
mix0.plotsubstancefit()

# ## Water

mix1 = SubstanceFit(T, P, {"H2O": 1.0})
mix1.plotsubstancefit()


# ## Combustion flue

def fluegases():
    """ Reference composition of flue gases used in model. """
    Y = {"CO2": 0.192, "O2": 0.016, "H2O": 0.076, "Ar": 0.012}
    Y["N2"] = 1.0 - sum(Y.values())
    return Y


mix2 = SubstanceFit(T, P, fluegases())
mix2.plotsubstancefit()


# ## Mixture law validation

class Mixture:
    """ Mixture quantities of substances. """
    def __init__(self, subs, T, Y):
        """ Create a mixture of substances. """
        if not np.allclose(sum(Y), 1.0):
            raise ValueError(f"Mass fractions must add up to 1.0: {Y}")

        self._sol = [s.to_solution(T=t) for s, t in zip(subs, T)]
        self._qty = [ct.Quantity(s, mass=y) for s, y in zip(self._sol, Y)]

        # Create a new object, Cantera is leaking (?!?!?)!
        self._mix = ct.Quantity(self._sol[0], mass=Y[0])

        # Quantity is not compatible with `sum` (?!?!?!)!
        for qk in self._qty[1:]:
            self._mix += qk

    @property
    def mixture(self):
        """ Access to internal mixture object. """
        return self._mix


# ### Case: isothermal mixture

tmix = 1000
mix = Mixture([mix0, mix1, mix2],
              [tmix, tmix, tmix],
              [ 0.2,  0.1,  0.7])

sol = mix.mixture
np.allclose(sol.cp_mass, sum([q.cp_mass * q.mass for q in mix._qty]))

(sol.viscosity, sum([q.viscosity * q.mass for q in mix._qty]))

(sol.thermal_conductivity, sum([q.thermal_conductivity * q.mass for q in mix._qty]))





# ## Create database

with open("mixtures.json", "w") as fp:
    json.dump({
        "mix0": mix0.coefficients,
        "mix1": mix1.coefficients,
        "mix2": mix2.coefficients
    }, fp, indent=4)


