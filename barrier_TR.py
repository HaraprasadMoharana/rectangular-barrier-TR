#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
barrier_TR.py
=============

Transmission and reflection probabilities of a one-dimensional rectangular
potential barrier, plotted against energy for a family of barrier-strength
parameters.

Run with:      python barrier_TR.py
Requires:      numpy, matplotlib
Author:        Haraprasad Moharana (@HaraprasadMoharana)
License:       MIT


-------------------------------------------------------------------------
1.  THE PHYSICAL PROBLEM
-------------------------------------------------------------------------

A particle of mass m and energy E > 0 travels from x = -inf toward a
barrier of height U0 and width L:

                         U(x)
                          ^
                          |        +-----------+
                     U0 --|        |           |
                          |        |           |
        --> incident      |        |  barrier  |
        <-- reflected     |        |           |     --> transmitted
                        0 +--------+-----------+-------------> x
                                   0           L

    U(x) = U0   for  0 <= x <= L,        U(x) = 0  elsewhere.

We solve the time-independent Schrodinger equation in the three regions and
match the wavefunction psi and its derivative dpsi/dx at x = 0 and x = L.
The result is a pair of probabilities, T and R, defined below.


-------------------------------------------------------------------------
2.  NOMENCLATURE  (every symbol used in this file)
-------------------------------------------------------------------------

  SYMBOL      CODE NAME     MEANING                                  UNITS
  ---------------------------------------------------------------------
  E           --            energy of the incident particle          J (or eV)
  U0          U0_eV         height of the potential barrier          J (or eV)
  L           L_nm          width (length) of the barrier            m (or nm)
  m           m_rel * M_E   particle mass; m_rel is the effective
                            mass in units of the free-electron mass  kg
  hbar        HBAR          reduced Planck constant, h / 2*pi        J s

  eps         eps           REDUCED ENERGY,  eps = E / U0.
                            Dimensionless. eps < 1 means the particle
                            is below the barrier top (classically
                            forbidden -> pure tunnelling); eps > 1
                            means it is above the barrier top.        --

  gamma       gamma         BARRIER STRENGTH PARAMETER (the quantity
                            this study sweeps),

                                gamma = (L / hbar) * sqrt(2 m U0)
                                gamma^2 = 2 m U0 L^2 / hbar^2

                            Dimensionless. It is the single number
                            that fixes the shape of the T(eps) curve:
                            height and width enter ONLY through this
                            combination. See section 3.                --

  k           --            wavenumber outside the barrier,
                            k = sqrt(2 m E) / hbar                     1/m
  kappa       x = kappa*L   DECAY CONSTANT inside the barrier for
                            eps < 1, kappa = sqrt(2 m (U0 - E)) / hbar.
                            Inside the barrier psi ~ exp(-kappa x),
                            so 1/kappa is the penetration depth.       1/m
  k'          y = k'*L      wavenumber inside the barrier for eps > 1,
                            k' = sqrt(2 m (E - U0)) / hbar             1/m

  T           T             TRANSMISSION PROBABILITY: the fraction of
                            incident probability current that emerges
                            at x > L. Equivalently, the probability
                            that one incident particle is found on the
                            far side. 0 <= T <= 1.                     --
  R           R             REFLECTION PROBABILITY: the fraction of
                            incident current that returns to x < 0.
                            0 <= R <= 1.                               --

  T + R = 1                 UNITARITY (conservation of probability
                            current). The barrier neither creates nor
                            absorbs particles. The script verifies
                            this numerically at every energy.

  n           n             resonance order, n = 1, 2, 3, ...          --

  D           D             shorthand used in the code for the
                            denominator factor 4*eps*|1-eps|; it has
                            no separate physical name.                 --


-------------------------------------------------------------------------
3.  WHAT "BARRIER STRENGTH" MEANS PHYSICALLY
-------------------------------------------------------------------------

Define the penetration (evanescent decay) length of a zero-energy particle
in the barrier:

        delta = hbar / sqrt(2 m U0)          [kappa at E = 0 is 1/delta]

Then

        gamma = L / delta

So gamma is simply the barrier width measured in units of the wavefunction's
natural decay length. It is a measure of OPACITY:

    gamma << 1   thin / weak barrier  -> the wavefunction barely decays
                 across it, T is large even well below the barrier top.
    gamma >> 1   thick / strong barrier -> the wavefunction is crushed
                 exponentially, T is exponentially small below the top and
                 the classical step behaviour is nearly recovered.

Two barriers with different U0 and L but the same gamma give IDENTICAL
T(eps) curves. That is why sweeping gamma (rather than U0 and L separately)
is the economical way to present this problem.


-------------------------------------------------------------------------
4.  CLOSED-FORM RESULT
-------------------------------------------------------------------------

  (a) eps < 1   -- TUNNELLING.  Inside the barrier the solution is a real
      exponential, and kappa*L = gamma*sqrt(1 - eps):

                                 sinh^2( gamma*sqrt(1-eps) )   -1
          T(eps) = [ 1  +  ------------------------------------- ]
                                     4 eps (1 - eps)

      Classically T would be exactly 0 here. It is not: this is quantum
      tunnelling. For a thick barrier the expression reduces to the
      familiar exponential law

          T  ~  16 eps (1 - eps) * exp( -2 gamma sqrt(1 - eps) ).

  (b) eps > 1   -- OVER-BARRIER SCATTERING.  Inside the barrier the solution
      oscillates, and k'*L = gamma*sqrt(eps - 1):

                                  sin^2( gamma*sqrt(eps-1) )    -1
          T(eps) = [ 1  +  ------------------------------------- ]
                                     4 eps (eps - 1)

      Classically T would be exactly 1 here. It is not: the particle can be
      reflected by a downward step. This is quantum reflection.

  (c) eps = 1   -- the two branches meet; the apparent 0/0 is removable:

          T(1) = 1 / (1 + gamma^2 / 4)

  (d) In all cases        R(eps) = 1 - T(eps).


-------------------------------------------------------------------------
5.  TRANSMISSION RESONANCES
-------------------------------------------------------------------------

Above the barrier T returns to exactly 1 whenever sin(k'L) = 0, i.e.

        k' L = n*pi        ->        eps_n = 1 + (n*pi / gamma)^2

Equivalently L = n * (lambda_inside / 2): an integer number of half
wavelengths fits inside the barrier, the waves reflected from the two edges
interfere destructively, and the barrier becomes perfectly transparent.
This is the exact quantum-mechanical analogue of an anti-reflection coating
or a Fabry-Perot cavity on resonance. The markers in panel (a) sit at these
energies. Larger gamma -> more resonances inside a given energy window.


-------------------------------------------------------------------------
6.  NUMERICAL NOTES
-------------------------------------------------------------------------

*  sinh(x) overflows 64-bit floats for x greater than about 710, which is
   reached at low eps once gamma is large. Above X_ASYMP the exact formula
   is replaced by its asymptotic form
        T -> 4 * [4 eps (1-eps)] * exp(-2x),        x = gamma*sqrt(1-eps)
   whose relative error is itself of order exp(-2x), hence negligible there.

*  R is computed from its own closed form, NOT as 1 - T. Near a resonance
   R collapses toward zero, and forming 1 - T there would lose all
   significant digits (catastrophic cancellation).

*  eps = 1 is handled as a special case within a narrow window TOL_EPS1.
"""

from __future__ import annotations

__author__  = "Haraprasad Moharana (@HaraprasadMoharana)"
__license__ = "MIT"

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ============================== CONFIG ====================================
# Edit these; everything else follows.

GAMMA_LIST = [1.0, 2.0, 5.0, 10.0, 20.0]   # barrier strength values to compare
EPS_MIN    = 1.0e-6                        # lowest reduced energy (avoid E = 0)
EPS_MAX    = 4.0                           # highest reduced energy (E = 4 U0)
N_EPS      = 6001                          # number of energy grid points

LOG_PANEL  = True        # also produce a log-scale T plot of the tunnelling region
OVERLAY    = True        # T and R overlaid on one axes, single gamma
THICKNESS  = True        # T and R against barrier thickness at fixed energy
MARK_RES   = True        # mark the over-barrier resonances on panel (a)

GAMMA_SHOW = 5.0         # gamma used by the overlaid figure
EPS_FIXED  = 1.5         # reduced energy used by the thickness figure (eps > 1)
SHOW_DEFS  = True        # print the symbol table when the script is run
SAVE       = True        # write PNG + PDF next to the script
OUTSTEM    = "barrier_TR"
DPI        = 600

C_T = "#1f6feb"          # colour for T in the overlaid figures
C_R = "#e8710a"          # colour for R in the overlaid figures

X_ASYMP    = 30.0        # kappa*L beyond which the asymptotic sinh form is used
TOL_EPS1   = 1.0e-13     # half-width of the eps = 1 special-case window

# Physical constants (SI), used only by gamma_from_physical()
HBAR = 1.054571817e-34      # reduced Planck constant           [J s]
M_E  = 9.1093837015e-31     # free-electron mass                [kg]
QE   = 1.602176634e-19      # elementary charge = J per eV      [C]


# ========================= DEFINITIONS PRINTOUT ===========================

DEFINITIONS = r"""
--------------------------------------------------------------------------
 SYMBOLS USED IN THIS SCRIPT
--------------------------------------------------------------------------
 E        incident particle energy
 U0       barrier height
 L        barrier width
 m        particle (or effective) mass
 hbar     reduced Planck constant

 eps      = E / U0          reduced energy          (dimensionless)
          eps < 1 : below the barrier top  -> tunnelling
          eps > 1 : above the barrier top  -> over-barrier scattering

 gamma    = (L/hbar)*sqrt(2 m U0)   barrier STRENGTH parameter
          = L / delta,  with delta = hbar/sqrt(2 m U0) the zero-energy
            penetration depth.  gamma is the barrier width expressed in
            units of the wavefunction decay length: a pure measure of
            how opaque the barrier is.  (dimensionless)

 kappa    = sqrt(2m(U0-E))/hbar   decay constant inside barrier (eps<1)
 k'       = sqrt(2m(E-U0))/hbar   wavenumber inside barrier     (eps>1)
          kappa*L = gamma*sqrt(1-eps)      k'*L = gamma*sqrt(eps-1)

 T        transmission probability: fraction of incident probability
          current that gets through to x > L
 R        reflection probability: fraction sent back to x < 0
 T + R    = 1  (unitarity / conservation of probability current)

 eps_n    = 1 + (n*pi/gamma)^2    over-barrier transmission resonances,
          where exactly n half-wavelengths fit inside the barrier and
          T = 1 (perfect transparency, Fabry-Perot condition)
--------------------------------------------------------------------------
"""


def print_definitions() -> None:
    """Print the symbol table to stdout (useful when sharing the script)."""
    print(DEFINITIONS)


# ========================= PUBLICATION STYLE ==============================

def set_style() -> None:
    """Nature / IEEE-flavoured rcParams: sans-serif, thin spines, inward ticks."""
    mpl.rcParams.update({
        "figure.dpi":          120,
        "savefig.dpi":         DPI,
        "savefig.bbox":        "tight",
        "savefig.pad_inches":  0.02,
        "font.family":         "sans-serif",
        "font.sans-serif":     ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":           9,
        "axes.labelsize":      10,
        "axes.titlesize":      10,
        "legend.fontsize":     8,
        "xtick.labelsize":     9,
        "ytick.labelsize":     9,
        "mathtext.fontset":    "dejavusans",
        "axes.linewidth":      0.8,
        "axes.labelpad":       3.0,
        "lines.linewidth":     1.4,
        "xtick.direction":     "in",
        "ytick.direction":     "in",
        "xtick.top":           True,
        "ytick.right":         True,
        "xtick.major.size":    3.5,
        "ytick.major.size":    3.5,
        "xtick.minor.size":    2.0,
        "ytick.minor.size":    2.0,
        "xtick.major.width":   0.8,
        "ytick.major.width":   0.8,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon":      False,
        "legend.handlelength": 1.6,
        "legend.labelspacing": 0.35,
    })


# ============================== PHYSICS ===================================

def transmission(eps, gamma: float) -> np.ndarray:
    """
    Transmission probability T of a rectangular barrier.

    T is the fraction of the incident probability current that emerges on the
    far side of the barrier, i.e. the probability that an incident particle is
    transmitted rather than reflected.

    Parameters
    ----------
    eps : array_like
        Reduced energy eps = E/U0 (dimensionless, >= 0).
    gamma : float
        Barrier strength parameter gamma = (L/hbar)*sqrt(2 m U0)
        (dimensionless). Equal to the barrier width in units of the
        zero-energy penetration depth.

    Returns
    -------
    ndarray
        T in [0, 1], same shape as eps.
    """
    eps = np.atleast_1d(np.asarray(eps, dtype=float))
    T   = np.empty_like(eps)

    below = eps < 1.0 - TOL_EPS1          # tunnelling regime
    above = eps > 1.0 + TOL_EPS1          # over-barrier regime
    at1   = ~(below | above)              # eps = 1 exactly

    # ---- eps < 1 : evanescent wave inside the barrier --------------------
    e = eps[below]
    x = gamma * np.sqrt(1.0 - e)          # x = kappa * L
    D = 4.0 * e * (1.0 - e)
    t = np.empty_like(e)
    ok = x < X_ASYMP
    t[ok]  = D[ok] / (D[ok] + np.sinh(x[ok]) ** 2)
    t[~ok] = 4.0 * D[~ok] * np.exp(-2.0 * x[~ok])    # sinh^2 x -> e^{2x}/4
    T[below] = t

    # ---- eps > 1 : propagating wave inside the barrier -------------------
    e = eps[above]
    y = gamma * np.sqrt(e - 1.0)          # y = k' * L
    D = 4.0 * e * (e - 1.0)
    T[above] = D / (D + np.sin(y) ** 2)

    # ---- eps = 1 : removable singularity ---------------------------------
    T[at1] = 1.0 / (1.0 + 0.25 * gamma ** 2)

    return T


def reflection(eps, gamma: float) -> np.ndarray:
    """
    Reflection probability R of a rectangular barrier.

    R is the fraction of incident probability current returned to the source
    side. Evaluated from its own closed form rather than as 1 - T, so that it
    stays accurate at the over-barrier resonances where R -> 0.

    Parameters and return shape are as in transmission().
    """
    eps = np.atleast_1d(np.asarray(eps, dtype=float))
    R   = np.empty_like(eps)

    below = eps < 1.0 - TOL_EPS1
    above = eps > 1.0 + TOL_EPS1
    at1   = ~(below | above)

    e = eps[below]
    x = gamma * np.sqrt(1.0 - e)
    D = 4.0 * e * (1.0 - e)
    r = np.empty_like(e)
    ok = x < X_ASYMP
    s2 = np.sinh(x[ok]) ** 2
    r[ok]  = s2 / (D[ok] + s2)
    r[~ok] = 1.0 - 4.0 * D[~ok] * np.exp(-2.0 * x[~ok])
    R[below] = r

    e = eps[above]
    y = gamma * np.sqrt(e - 1.0)
    D = 4.0 * e * (e - 1.0)
    s2 = np.sin(y) ** 2
    R[above] = s2 / (D + s2)

    R[at1] = 1.0 - 1.0 / (1.0 + 0.25 * gamma ** 2)

    return R


def resonance_energies(gamma: float, eps_max: float) -> np.ndarray:
    """
    Over-barrier transmission resonances eps_n = 1 + (n*pi/gamma)^2, n >= 1,
    lying below eps_max. At these energies an integer number of half
    wavelengths fits inside the barrier and T = 1 exactly.
    """
    n_max = int(np.floor(gamma * np.sqrt(max(eps_max - 1.0, 0.0)) / np.pi))
    n = np.arange(1, n_max + 1)
    return 1.0 + (n * np.pi / gamma) ** 2


def gamma_from_physical(U0_eV: float, L_nm: float, m_rel: float = 1.0) -> float:
    """
    Convert real barrier parameters into the dimensionless strength

        gamma = (L / hbar) * sqrt(2 m U0)

    Parameters
    ----------
    U0_eV : barrier height in electron-volts
    L_nm  : barrier width in nanometres
    m_rel : effective mass in units of the free-electron mass m_e
            (use 1.0 for a free electron)

    Returns
    -------
    float : gamma (dimensionless)
    """
    L  = L_nm * 1.0e-9
    U0 = U0_eV * QE
    m  = m_rel * M_E
    return L * np.sqrt(2.0 * m * U0) / HBAR


def penetration_depth_nm(U0_eV: float, m_rel: float = 1.0) -> float:
    """
    Zero-energy penetration depth delta = hbar / sqrt(2 m U0), in nanometres.
    gamma is just L / delta.
    """
    m = m_rel * M_E
    return 1.0e9 * HBAR / np.sqrt(2.0 * m * U0_eV * QE)


# ============================== PLOTTING ==================================

def make_main_figure(eps: np.ndarray, gammas, mark_res: bool = True):
    """Side-by-side T(eps) and R(eps), one curve per barrier strength gamma."""
    colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(gammas)))

    fig, (axT, axR) = plt.subplots(1, 2, figsize=(7.2, 3.1), sharex=True)

    for g, c in zip(gammas, colors):
        lab = rf"$\gamma = {g:g}$"
        axT.plot(eps, transmission(eps, g), color=c, label=lab)
        axR.plot(eps, reflection(eps, g),  color=c, label=lab)

    if mark_res:
        for g, c in zip(gammas, colors):
            for er in resonance_energies(g, eps.max()):
                axT.plot(er, 1.0, marker="o", ms=2.6, mfc="none",
                         mec=c, mew=0.7, zorder=5, clip_on=False)

    for ax in (axT, axR):
        ax.axvline(1.0, color="0.55", lw=0.7, ls=(0, (4, 3)), zorder=0)
        ax.set_xlim(0.0, eps.max())
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel(r"Reduced energy  $\varepsilon = E/U_0$")

    axT.set_ylabel(r"Transmission  $T$")
    axR.set_ylabel(r"Reflection  $R$")
    axT.text(0.03, 0.94, "(a)", transform=axT.transAxes, va="top", fontweight="bold")
    axR.text(0.03, 0.94, "(b)", transform=axR.transAxes, va="top", fontweight="bold")

    axT.legend(loc="lower right", title=r"$\gamma=\frac{L}{\hbar}\sqrt{2mU_0}$",
               title_fontsize=8)

    fig.tight_layout(w_pad=1.8)
    return fig


def make_log_figure(gammas):
    """
    Tunnelling region only, on a logarithmic T axis. The near-straight
    families expose the exp(-2 gamma sqrt(1-eps)) law.
    """
    eps = np.linspace(EPS_MIN, 1.0 - 1e-9, 4000)
    colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(gammas)))

    fig, ax = plt.subplots(figsize=(3.6, 3.1))
    for g, c in zip(gammas, colors):
        ax.semilogy(eps, transmission(eps, g), color=c, label=rf"$\gamma = {g:g}$")

    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"Reduced energy  $\varepsilon = E/U_0$")
    ax.set_ylabel(r"Transmission  $T$")
    ax.legend(loc="lower right")
    fig.tight_layout()
    return fig


def transmission_vs_thickness(u, eps: float) -> np.ndarray:
    """
    Transmission as a function of barrier THICKNESS at fixed reduced energy.

    This is the perpendicular cut through the same closed form used by
    transmission(): there the barrier is fixed and the energy sweeps, here the
    energy is fixed and the barrier grows. The two run in opposite directions --
    a thicker barrier transmits less, a faster particle transmits more -- so the
    curves look inverted even though the physics is identical.

    Parameters
    ----------
    u : array_like
        Dimensionless barrier thickness. For eps > 1 this is k'L; for eps < 1
        it is kappa*L. In both cases u = gamma*sqrt(|1 - eps|).
    eps : float
        Fixed reduced energy E/U0.

    Returns
    -------
    ndarray
        T in [0, 1], same shape as u. At u = 0 (no barrier) T = 1 exactly.
    """
    u = np.atleast_1d(np.asarray(u, dtype=float))
    if eps > 1.0:
        D = 4.0 * eps * (eps - 1.0)
        return D / (D + np.sin(u) ** 2)
    if eps < 1.0:
        D = 4.0 * eps * (1.0 - eps)
        return D / (D + np.sinh(np.clip(u, 0.0, X_ASYMP)) ** 2)
    return 1.0 / (1.0 + 0.25 * u ** 2)          # eps = 1 limit


def _strip(ax) -> None:
    """Drop the top and right spines and their ticks."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(which="both", top=False, right=False)


def make_overlaid_figure(gamma: float = GAMMA_SHOW):
    """
    T and R for a single gamma, overlaid on one axes against reduced energy.

    The two-panel figure separates T and R; this one puts them together so the
    crossing at T = R = 1/2 and the unitarity constraint T + R = 1 are visible
    directly. The shaded strip is the classically forbidden region eps < 1.
    """
    eps = np.linspace(EPS_MIN, EPS_MAX, N_EPS)
    T, R = transmission(eps, gamma), reflection(eps, gamma)

    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.axvspan(0.0, 1.0, color="0.94", zorder=0, lw=0)
    ax.plot(eps, T, color=C_T, lw=1.9, zorder=3)
    ax.plot(eps, R, color=C_R, lw=1.9, zorder=3)

    for er in resonance_energies(gamma, EPS_MAX):
        ax.plot(er, 1.0, marker="o", ms=4, mfc="white", mec=C_T, mew=1.1, zorder=5)

    ax.axvline(1.0, color="0.6", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(0.875 * EPS_MAX, 0.84, "$T$", color=C_T, fontsize=13,
            fontweight="bold", ha="center")
    ax.text(0.875 * EPS_MAX, 0.13, "$R$", color=C_R, fontsize=13,
            fontweight="bold", ha="center")
    ax.text(0.5, 1.05, "tunnelling", ha="center", fontsize=8.5, color="0.4")
    ax.text(0.5 * (1.0 + EPS_MAX), 1.05, "over-barrier", ha="center",
            fontsize=8.5, color="0.4")

    ax.set_xlim(0.0, EPS_MAX)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"Reduced energy  $\varepsilon = E/U_0$")
    ax.set_ylabel("Probability")
    ax.set_title(rf"$\gamma = {gamma:g}$", fontsize=10, pad=18)
    _strip(ax)
    fig.tight_layout()
    return fig


def make_thickness_figure(eps0: float = EPS_FIXED, n_periods: int = 4):
    """
    T and R against barrier thickness at fixed energy, for eps0 > 1.

    T starts at 1 for a vanishing barrier and dips, returning to 1 at every
    k'L = n*pi where an integer number of half wavelengths fits inside. This is
    the Fabry-Perot / anti-reflection condition seen along the thickness axis
    instead of the energy axis.
    """
    if eps0 <= 1.0:
        raise ValueError("make_thickness_figure expects eps0 > 1 (over-barrier)")

    u = np.linspace(0.0, n_periods * np.pi, 4000)
    T = transmission_vs_thickness(u, eps0)
    R = 1.0 - T

    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.plot(u, T, color=C_T, lw=1.9)
    ax.plot(u, R, color=C_R, lw=1.9)
    for n in range(1, n_periods + 1):
        ax.plot(n * np.pi, 1.0, marker="o", ms=4, mfc="white",
                mec=C_T, mew=1.1, zorder=5)

    ax.text(np.pi / 2, 0.64, "$T$", color=C_T, fontsize=13,
            fontweight="bold", ha="center")
    ax.text(np.pi / 2, 0.36, "$R$", color=C_R, fontsize=13,
            fontweight="bold", ha="center")

    ax.set_xlim(0.0, n_periods * np.pi)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([n * np.pi for n in range(n_periods + 1)])
    ax.set_xticklabels(["0", r"$\pi$"] +
                       [rf"${n}\pi$" for n in range(2, n_periods + 1)])
    ax.set_xlabel(r"Barrier thickness  $k'L$")
    ax.set_ylabel("Probability")
    ax.set_title(rf"$\varepsilon = {eps0:g}$   (over-barrier, $E > U_0$)",
                 fontsize=10, pad=18)
    _strip(ax)
    fig.tight_layout()
    return fig


# =============================== MAIN =====================================

def main() -> None:
    if SHOW_DEFS:
        print_definitions()

    set_style()
    eps = np.linspace(EPS_MIN, EPS_MAX, N_EPS)

    # Unitarity check: T + R - 1 should sit at machine precision everywhere.
    worst = 0.0
    for g in GAMMA_LIST:
        dev = np.max(np.abs(transmission(eps, g) + reflection(eps, g) - 1.0))
        worst = max(worst, dev)
    print(f"unitarity check:  max |T + R - 1|  =  {worst:.3e}")

    fig1 = make_main_figure(eps, GAMMA_LIST, mark_res=MARK_RES)
    if SAVE:
        fig1.savefig(f"{OUTSTEM}_vs_eps.png")
        fig1.savefig(f"{OUTSTEM}_vs_eps.pdf")

    if LOG_PANEL:
        fig2 = make_log_figure(GAMMA_LIST)
        if SAVE:
            fig2.savefig(f"{OUTSTEM}_log.png")
            fig2.savefig(f"{OUTSTEM}_log.pdf")

    if OVERLAY:
        fig3 = make_overlaid_figure(GAMMA_SHOW)
        if SAVE:
            fig3.savefig(f"{OUTSTEM}_overlaid.png")
            fig3.savefig(f"{OUTSTEM}_overlaid.pdf")

    if THICKNESS:
        fig4 = make_thickness_figure(EPS_FIXED)
        if SAVE:
            fig4.savefig(f"{OUTSTEM}_thickness.png")
            fig4.savefig(f"{OUTSTEM}_thickness.pdf")

    plt.show()


if __name__ == "__main__":
    main()
