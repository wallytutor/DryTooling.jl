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

        # Store internally the solution array for conversion later.
        # Since all elements have the same `mw`, keep first only.
        self._sol = sol
        self._mw = sol.mean_molecular_weight[1]
        self._Y = Y

        def fitrefs(p, d):
            return np.polyfit(sol.T, self._refs[p], deg=d)

        self._refs = {p: getattr(sol, p) for p, d in degs.items()}
        self._fits = {p: fitrefs(p, d)   for p, d in degs.items()}

        # Store coefficients with enough infor for reproducible builds.
        # Mean molecular weights are required for PFR mixing laws.
        self._coef = {
            "mech": mech,
            "Y": self._Y,
            "W": self._mw,
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
        xtick = np.arange(T[0], T[-1]+1, dT)

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

    def specific_heat(self, T):
        """ Evaluate substance specific heat through polynomial. """
        return np.polyval(self._fits["cp_mass"], T)

    def to_solution(self, T=1000.0, P=P):
        """ Get standard Cantera solution from mixture. """
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

# The following formulas from Kee et al. (2017) are proposed for what follows:
#
# - Viscosity (Wilke-Bird):
#
# $$
# \mu = \sum_{k=1}^{K}\dfrac{X_k\mu_k}{\sum_{j=1}^{K}X_j\Phi_{kj}}
# $$
#
# $$
# \Phi_{kj}=\frac{1}{\sqrt{8}}
# \left(1+\frac{W_k}{W_j}\right)^{-1/2}
# \left[1+\left(\frac{\mu_k}{\mu_j}\right)^{1/2}\left(\frac{W_j}{W_k}\right)^{1/4}\right]^2
# $$
#
# - Thermal conductivity (Mathur et al.):
#
# $$
# \lambda = \frac{1}{2}\left(\sum_{k=1}^{K}X_k\lambda_k + \dfrac{1}{\sum_{k=1}^{K}X_k/\lambda_k}\right)
# $$

class Mixture:
    """ Mixture quantities of substances. """
    def __init__(self, subs, T, Y):
        """ Create a mixture of substances. """
        if not np.allclose(sum(Y), 1.0):
            raise ValueError(f"Mass fractions must add up to 1.0: {Y}")

        self._sub = subs
        self._T = T

        self._sol = [s.to_solution(T=t) for s, t in zip(subs, T)]
        self._qty = [ct.Quantity(s, mass=y) for s, y in zip(self._sol, Y)]

        self._mw = np.array([q.mean_molecular_weight for q in self._qty])
        self._mu = np.array([q.viscosity for q in self._qty])
        self._kg = np.array([q.thermal_conductivity for q in self._qty])

        # Create a new object, Cantera is leaking (?!?!?)!
        self._mix = ct.Quantity(self._sol[0], mass=Y[0])

        # Quantity is not compatible with `sum` (?!?!?!)!
        for qk in self._qty[1:]:
            self._mix += qk

        self._Y = np.array(Y)
        self._X = self._Y * self._mw / self._mix.mean_molecular_weight

    @property
    def Y(self):
        """ Access mass fractions of solutions. """
        return self._Y

    @property
    def X(self):
        """ Access mole fractions of solutions. """
        return self._X

    @property
    def mean_molecular_weight(self):
        """ Access mixture mean molecular mass. """
        return self._mix.mean_molecular_weight

    @property
    def mixture(self):
        """ Access to internal mixture object. """
        return self._mix

    def mass_weighted_specific_heat(self, usepoly=False):
        """ Mass weighted specific heat model. """
        # For isothermal mixtures this would be enough.
        # return sum([q.cp_mass * q.mass for q in self._qty])
        K = len(mix._Y)

        if not usepoly:
            cps = [self._qty[k].cp_mass for k in range(K)]
        else:
            cps = [self._sub[k].specific_heat(self._T[k])  for k in range(K)]

        den = [self._Y[k] * cps[k]  for k in range(K)]
        num = [den[k] * self._qty[k].T for k in range(K)]
        Tm = sum(num) / sum(den)

        cp = 0
        for k, s in enumerate(self._sol):
            if not usepoly:
                Told = s.T
                s.TP = Tm, None
                cp += s.cp_mass * self._Y[k]
                s.TP = Told, None
            else:
                f = self._sub[k].specific_heat
                cp += f(Tm) * self._Y[k]

        return cp

    def mass_weighted_viscosity(self):
        """ Mass weighted viscosity model. """
        return sum([q.viscosity * q.mass for q in self._qty])

    def mass_weighted_thermal_conductivity(self):
        """ Mass weighted thermal conductivity model. """
        return sum([q.thermal_conductivity * q.mass for q in self._qty])

    def mole_weighted_viscosity(self):
        """ Mole weighted viscosity model. """
        return sum([q.viscosity * self._X[k]
                    for k, q in enumerate(self._qty)])

    def mole_weighted_thermal_conductivity(self):
        """ Mole weighted thermal conductivity model. """
        return sum([q.thermal_conductivity * self._X[k]
                    for k, q in enumerate(self._qty)])

    def wilke_bird_viscosity(self):
        """ Wilke-Bird viscosity model. """
        X, W, mu = self._X, self._mw, self._mu

        def phi_kj(muk, muj, wk, wj):
            s = (1 + (muk/muj)**(1/2) * (wj/wk)**(1/4))**2
            return s / (8 * (1 + wk/wj))**(1/2)

        m = 0
        for k in range(len(X)):
            d = sum(X[j] * phi_kj(mu[k], mu[j], W[k], W[j])
                    for j in range(len(X)))
            m += (X[k] * mu[k] / d)

        return m

    def mathur_thermal_condictivity(self):
        """ Mathur et al. thermal conductivity model. """
        X, k = self._X, self._kg
        return (1/2)*((X @ k) + 1 / sum(X / k))


# ### Case: isothermal mixture

# +
tmix = [500]*3
comp = [0.1, 0.1, 0.8]
mix = Mixture([mix0, mix1, mix2], tmix, comp)

print(mix.mixture.cp_mass,
      mix.mass_weighted_specific_heat())

print(mix.mixture.viscosity,
      mix.mass_weighted_viscosity(),
      mix.mole_weighted_viscosity(),
      mix.wilke_bird_viscosity())

print(mix.mixture.thermal_conductivity,
      mix.mass_weighted_thermal_conductivity(),
      mix.mole_weighted_thermal_conductivity(),
      mix.mathur_thermal_condictivity())
# -

# Now we test it over a broad range of temperatures of interest.

temperatures = np.arange(300, 2001, 100)
mixtures = [Mixture([mix0, mix1, mix2], [tk]*3, comp)
            for tk in temperatures]

# We recover arrays of properties in a organized format.

# +
mus = np.array([
        [mix.mixture.viscosity,
         mix.mass_weighted_viscosity(),
         mix.mole_weighted_viscosity(),
         mix.wilke_bird_viscosity()]
         for mix in mixtures
        ])

kgs = np.array([
        [mix.mixture.thermal_conductivity,
         mix.mass_weighted_thermal_conductivity(),
         mix.mole_weighted_thermal_conductivity(),
         mix.mathur_thermal_condictivity()]
         for mix in mixtures
        ])


# -

# Finally we compute the mean relative squared error of the models...

# +
def mrse(a, b):
    """ Metric for model performance comparison. """
    return sum(((a - b) / a)**2)


mus_err = (mrse(mus[:, 0], mus[:, 1]),
           mrse(mus[:, 0], mus[:, 2]),
           mrse(mus[:, 0], mus[:, 3]))

kgs_err = (mrse(kgs[:, 0], kgs[:, 1]),
           mrse(kgs[:, 0], kgs[:, 2]),
           mrse(kgs[:, 0], kgs[:, 3]))
# -

# ...and display the results for taking a conclusion.

# +
plt.close("all")
plt.style.use("seaborn-white")
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax[0].plot(temperatures, mus[:, 0], "ko", label="Reference")
ax[1].plot(temperatures, kgs[:, 0], "ko", label="Reference")

ax[0].plot(temperatures, mus[:, 1], "r", label=f"Mass weighted {mus_err[0]:.3e}")
ax[1].plot(temperatures, kgs[:, 1], "r", label=f"Mass weighted  {kgs_err[0]:.3e}")

ax[0].plot(temperatures, mus[:, 2], "b", label=f"Mole weighted {mus_err[1]:.3e}")
ax[1].plot(temperatures, kgs[:, 2], "b", label=f"Mole weighted {kgs_err[1]:.3e}")

ax[0].plot(temperatures, mus[:, 3], "g", label=f"Model {mus_err[2]:.3e}")
ax[1].plot(temperatures, kgs[:, 3], "g", label=f"Model {kgs_err[2]:.3e}")

ax[0].grid(linestyle=":")
ax[1].grid(linestyle=":")

ax[0].set_ylabel("Viscosity [Pa.s]")
ax[1].set_ylabel("Conductivity [W/(m.K)]")

ax[0].legend(loc=4)
ax[1].legend(loc=4)

ax[0].set_xlim(300, 2000)
ax[0].set_xticks(range(300, 2001, 100))

fig.tight_layout()
# -

# Because of good performance in both cases, a mass weighted model is recommended for use with PFR's under a simplified framework. Furthermore, this is the least expensive model to be computed.

# ### Case: non-isothermal mixture

# +
Ya = 0.1

tmix = [500, 1000]
comp = [Ya, 1-Ya]
mixs = [mix0, mix2]
mix = Mixture(mixs, tmix, comp)

print(mix.mixture.cp_mass,
      mix.mass_weighted_specific_heat(usepoly=True))

# +
# TODO scan composition.
# -

# ## Create database

with open("mixtures.json", "w") as fp:
    json.dump({
        "mix0": mix0.coefficients,
        "mix1": mix1.coefficients,
        "mix2": mix2.coefficients
    }, fp, indent=4)
