#!/bin/env python2.7
# coding: utf-8

import os.path
from datetime import datetime

import click
import numpy as np
from scipy.optimize import curve_fit

_e = 1.60217733e-19
kJ = 1000.0 / _e


class EOS(object):
    """ Generic function for eos fitting. You need to subclass this and redefine
    fitting functions. Dump_fit will give you a range of points fitted using
    results from your data
    """

    EOS_NAME = "Generic EOS"

    def __init__(self, volumes, energies):
        self.v = np.array(volumes)
        self.e = np.array(energies)

        self.v0 = None

    @classmethod
    def parabola(cls, x, a, b, c):
        """
        parabola polynomial function

        this function is used to fit the data to get good guesses for
        the equation of state fits

        a 4th order polynomial fit to get good guesses for
        was not a good idea because for noisy data the fit is too wiggly
        2nd order seems to be sufficient, and guarentees a single minimum"""

        return a + b * x + c * x ** 2

    @classmethod
    def energy_fit_function(cls, V, E0, B0, BP, V0):
        pass

    @classmethod
    def pressure_fit_function(cls, V, E0, B0, BP, V0):
        pass

    def fit(self):
        """Calculate volume, energy, and bulk modulus.
        old ASE2 implementation

        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the ASE units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          v0, e0, B, BP = eos.fit()
          print(B / kJ * 1.0e24, 'GPa')

        """

        p0 = [min(self.e), 1, 1]
        popt, pcov = curve_fit(self.parabola, self.v, self.e, p0)

        parabola_parameters = popt
        # Make sure the minimum is bracketed by the volumes
        # (this if for the solver)
        minvol = min(self.v)
        maxvol = max(self.v)

        # the minimum of the parabola is at dE/dV = 0, or 2*c V +b =0
        c = parabola_parameters[2]
        b = parabola_parameters[1]
        a = parabola_parameters[0]
        parabola_vmin = -b / 2 / c

        if not (minvol < parabola_vmin and parabola_vmin < maxvol):
            print(
                "Warning the minimum volume of a fitted parabola is not in your volumes. You may not have a minimum in your dataset"
            )

        # evaluate the parabola at the minimum to estimate the groundstate energy
        E0 = self.parabola(parabola_vmin, a, b, c)
        # estimate the bulk modulus from Vo*E''.  E'' = 2*c
        B0 = 2 * c * parabola_vmin

        BP = 4  # just a guess

        initial_guess = [E0, B0, BP, parabola_vmin]

        # now fit the equation of state
        p0 = initial_guess
        popt, pcov = curve_fit(self.energy_fit_function, self.v, self.e, p0)

        self.eos_parameters = popt

        self.v0 = self.eos_parameters[3]
        self.e0 = self.eos_parameters[0]
        self.B = self.eos_parameters[1]
        self.BP = self.eos_parameters[2]

        return self.v0, self.e0, self.B, self.BP

    def dump_fit(self, filename):
        """ Output table of volumes into file, energies and pressures to stdout.
        Uses the volumes provided for fitting +- 10%

        """
        if self.v0 is None:
            self.fit()

        vol_list = np.linspace(min(self.v) * 0.95, max(self.v) * 1.1)
        vol_list = np.sort(np.append(vol_list, self.v))

        def pressure(v):
            return self.pressure_fit_function(
                v,
                self.eos_parameters[0],
                self.eos_parameters[1],
                self.eos_parameters[2],
                self.eos_parameters[3],
            )

        def energy(v):
            return self.energy_fit_function(
                v,
                self.eos_parameters[0],
                self.eos_parameters[1],
                self.eos_parameters[2],
                self.eos_parameters[3],
            )

        with open(filename, "w") as fp:
            now_str = datetime.now().strftime("%Y-%b-%d %H:%M")
            print("# %s EOS fit %s" % (self.EOS_NAME, now_str), file=fp)
            print("# V[Ang^3] Energy[eV] P[GPa]", file=fp)
            for v, e, p in zip(vol_list, energy(vol_list), pressure(vol_list)):
                print(v, e, p / kJ * 1e24, file=fp)


class BirchMurnaghanEOS(EOS):
    """ Fit equation of state for bulk systems.
    Based on EquationOfStateASE2, but adds BP and fixed to Birch-Murnaghan

       birchmurnaghan
           as seen in PRB 70, 224107
           Original paper by Birch: DOI: 10.1103/PhysRev.71.809

    Use::

       eos = BirchMurnaghanEOS(volumes, energies)
       v0, e0, B, BP = eos.fit()

    """

    EOS_NAME = "Birch-Murnaghan"

    @classmethod
    def energy_fit_function(cls, V, E0, B0, BP, V0):
        """
        BirchMurnaghan equation from PRB 70, 224107
        Those fuckers used wrong ETA there. Bitches
        """

        eta = (V0 / V) ** (2.0 / 3.0)
        E = E0 + 9.0 * B0 * V0 / 16.0 * (eta - 1) ** 2 * (
            6 + BP * (eta - 1.0) - 4.0 * eta
        )
        return E

    @classmethod
    def pressure_fit_function(cls, V, E0, B0, BP, V0):
        "B-M pressure as a function of Volume using supplied parameters"

        eta = (V0 / V) ** (1.0 / 3.0)
        P = (
            3.0
            * B0
            / 2.0
            * (eta ** 7 - eta ** 5)
            * (1.0 + 3.0 * (BP - 4.0) / 4.0 * (eta ** 2 - 1))
        )
        return P


class VinetEOS(EOS):
    """ Fit equation of state for bulk systems.
    Based on EquationOfStateASE2, but adds BP and fixed to vinet

       vinet
           PRB 70, 224107

    Use::

       eos = VinetEOS(volumes, energies)
       v0, e0, B, BP = eos.fit()

    """

    EOS_NAME = "Vinet"

    @classmethod
    def energy_fit_function(cls, V, E0, B0, BP, V0):
        "Vinet equation from PRB 70, 224107"

        eta = (V / V0) ** (1.0 / 3.0)

        E = E0 + 2.0 * B0 * V0 / (BP - 1.0) ** 2 * (
            2.0
            - (5.0 + 3.0 * BP * (eta - 1.0) - 3.0 * eta)
            * np.exp(-3.0 * (BP - 1.0) * (eta - 1.0) / 2.0)
        )
        return E

    @classmethod
    def pressure_fit_function(cls, V, E0, B0, BP, V0):
        "Vinet pressure as a function of Volume using supplied parameters"

        eta = (V / V0) ** (1.0 / 3.0)
        P = (
            3.0
            * B0
            * (1.0 - eta)
            / (eta ** 2)
            * np.exp(3.0 / 2.0 * (BP - 1.0) * (1.0 - eta))
        )
        return P


@click.command()
@click.option("--datafile", prompt="Please provide datafile name", help="datafile name")
def main(datafile):
    """ Reads lines from `datafile`, lines expected to be in format

    volume/cell energy/cell

    You need at least a few points around minimum to make the fit trustworthy

    """

    d = []
    with open(datafile, "r") as fp:
        for line in fp:
            if not line.startswith("#"):
                s = line.split()
                d += [(float(s[0]), float(s[1]))]

    eos = VinetEOS([v for v, _ in d], [e for _, e in d])
    V0, E0, B, BP = eos.fit()
    B = B / kJ * 1.0e24
    print("vinet equation of state fit results:")
    print("  V0 = %s Ang^3" % V0)
    print("  E0 = %s eV" % E0)
    print("  B = %s GPa" % B)
    print("  BP = %s" % BP)

    plot_fname = "%s.vinet_plot" % os.path.splitext(datafile)[0]
    eos.dump_fit(plot_fname)
    print("plot data saved as '%s'" % plot_fname)

    print("-".center(60, "-"))

    eos = BirchMurnaghanEOS([v for v, _ in d], [e for _, e in d])
    V0, E0, B, BP = eos.fit()
    B = B / kJ * 1.0e24
    print("Birch-Murnaghan equation of state fit results:")
    print("  V0 = %s Ang^3" % V0)
    print("  E0 = %s eV" % E0)
    print("  B = %s GPa" % B)
    print("  BP = %s" % BP)

    plot_fname = "%s.bm_plot" % os.path.splitext(datafile)[0]
    eos.dump_fit(plot_fname)
    print("plot data saved as '%s'" % plot_fname)


if __name__ == "__main__":
    main()
