#!/usr/bin/python
#=============================================================================================
# Constants.py
#
# Mathematical and physical constants.
#
# Written 2007-05-17
# John D. Chodera <jchodera@gmail.com>
# Department of Chemistry, Stanford University
# Pande group
#
#Inspired by Units.py from MMTK, by Konrad Hinsen.
#=============================================================================================
#CVS VERSION CONTROL INFORMATION
__version__ = "$Revision: $"
# $Source: $
#=============================================================================================
# TODO
# - Provide complete names and aliases in addition to symbolic names.
# - Is there a way to ensure each constant is documented by help() or pydoc?
# - Complete documentation.
# - Refer to NIST for updated constants?
# - Make sure only derives from Units -- eliminate cross-dependencies.
# - Give full names in addition to common abbreviations for all constants.
#=============================================================================================
"""
Mathematical and physical constants.

Mathematical constants are imported from the "math" package, and generally provided for convenience.

All physical constants come from the CRC Handbook of Chemistry and Physics, 59th Edition (1978-1979),
"Recommended consistentt values of the fundamental physical constants", page F-250.

Provided constants:

c (speed of light),
Na (Avogadro's number),
h = (Planck constant),
hbar = (Planck constant divided by 2*Pi),
k_B = (Boltzmann constant),
eps0 = (permittivity of vacuum),
me = (electron mass)
"""

#=============================================================================================
# MATHEMATICAL CONSTANTS
#=============================================================================================

# pi
from math import pi

#=============================================================================================
# PHYSICAL CONSTANTS
#=============================================================================================
import Units

# Avogadro number
from Units import N_A # Avogadro's number

# permeability of vacuum
mu_0 = 4.0 * pi * 1.0e-7 * Units.henry / Units.meter

# speed of light in vacuum
c = 2.99792458e8 * Units.m / Units.s**2

# permittivity of vacuum
epsilon_0 = 1.0 / (mu_0 / c**2)

# elementary charge
e = 1.6021892e-19 * Units.Coulomb

# Planck constant
h = 6.626176e-34 * Units.J * Units.s
hbar = h / (2.0 * pi)

# fine structure constant
alpha = (mu_0 * c * e**2) / (2.0 * h)

# electron rest mass
m_e = 5.4858026e-4 * Units.amu

# muon rest mass
m_mu = 0.11342920 * Units.amu

# proton rest mass
m_p = 1.007276470 * Units.amu

# neutron rest mass
m_n = 1.008665012 * Units.amu

# Faraday constant
Faraday_constant = N_A * e

# Rydberg constant
R_infinity = 1.097373177e7 / Units.m

# Bohr radius
alpha_0 = alpha / (4.0 * pi * R_infinity)
Bohr = alpha_0 # alias for Bohr radius

# Bohr magneton
mu_B = e * hbar / (2.0 * m_e)

# nuclear magneton
mu_N = e * hbar / (2.0 * m_p)

# electron magnetic moment
mu_e = 9.284832e-24 * Units.joule / Units.tesla

# proton magnetic moment
mu_p = 1.4106171e-26 * Units.joule / Units.tesla

# muon magnetic moment
mu_mu = 4.490474e-26 * Units.joule / Units.tesla

# proton gyromagnetic ratio
gamma_p = 2.6751987e8 / Units.second / Units.tesla

# molar gas constant
R = 8.31441 * Units.J / Units.mol / Units.K

# molar volume, ideal gas (T_0 = 273.15 K, p_0 = 1 atm)
V_m = R * 273.15 * Units.K / (1.0 * Units.atm)

# Boltzmann constant
kB = R / N_A
k_B = kB # alias

# gravitational constant
G = 6.6720e-11 * Units.N * Units.m**2 / Units.kg**2

# Hartree
Hartree = hbar**2 / (m_e * Bohr**2) 
