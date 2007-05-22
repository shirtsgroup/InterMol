#!/usr/bin/python
#=============================================================================================
# Units.py
#
# Utility class to facilitate working with unit-attached numericial quantities.

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
# - Fix dependencies for constants N_A and pi.
# - Update documentation.
# - Refer to NIST
# - Give full names in addition to common abbreviations for all units.
# - Split into separate scopes for SI, MKS, CGS, English, Common?
# - Double-check units.
#=============================================================================================
"""
Utility class to facilitate working with unit-attached numerical quantities.

USE

To associate units with a quantity, simply mutiply by the appropriate units:

  import Units
  # Define a box edge length in Angstroms.
  length = 10.5 * Units.A
  # Define a velocity.
  velocity = 4.5 * Units.m / Units.s # in m/s
  # Define a force.
  force = 3.5 * Units.N # in Newtons

To extract quantity in some desired set of units, simply divide by the desired units:

  import Units
  # Extract box length in nanometers.
  length_in_nm = length / Units.nm
  # Extract the velocity in cm/s.
  velocity_in_cm_per_s = velocity / (Units.cm/Units.s)
  # Extract the force in pN/nm.
  force_in_pN_per_nm = force / (Units.pN/Units.nm)
  
WARNINGS

Once units are attached to a quantity, failure to divide by the units desired will likely result in erroneous behavior.
Quantities that have units attached are represented in an internal standard set of units, which may not be the units desired
and may even change from version to version.

  # Don't do this.
  D = 4.82e-5 Units.cm**2 / Units.s # diffusion coefficient of oxygen in water at 60 C
  t = 3.56 * Units.us # observation time  
  mean_squared_displacement = 6.0 * D * t**2
  print "Mean-squared displacement is %f." % mean_squared_displacement # forgot to convert back to sensible units!!!

Mixing unit-attached quantities and non-unit-attached quantities (that are not dimensionless) will get you into trouble.

  # Don't do this either.
  D = 4.82e-5 Units.cm**2 / Units.s # diffusion coefficient of oxygen in water at 60 C
  t = 3.56 # forgot to attache units!!!
  mean_squared_displacement = 6.0 * D * t**2 # this won't make any sense anymore

No sophisticated unit-tracking is done, so if you try to convert into something that doesn't make sense dimensionally, no error will be raised.

  # Especially don't do this.
  energy = 1.34 * Units.joules
  my_energy = energy / (Units.kg * Units.m**3 / Units.s) # denominator not in units of energy!

Don't use "from Units import *" to import the Units class -- this will just barf all over your local namesace in a most unpleasant manner.

  # Please don't do this.
  from Units import *
  # Always do this instead
  import Units

EXAMPLES

An example of computing the potential energy of a person standing on a wall:

  # Compute the potential energy due to gravity.
  import Units
  m = 42.0 * Units.kg # mass of object
  g = 9.8 * Units.m / Units.s**2 # acceleration due to gravity at sea level
  h = 3.0 * Units.m # height of object
  V = m * g * h # compute potential energy
  print "Potential energy is %f Joules" % V / Units.J # print energy in desired units

An simple illustration of velocity Verlet:

  # An example of velocity Verlet.
  dt = 1.0 * Units.fs # timestep size
  tmax = 1.0 * Units.ns # quantity of time to integrate
  t = 0.0 * Units.ps # current time counter
  m = 12.0 * Units.amu # mass of particle
  # Initialize position and velocities.
  (x_new, v_new) = (1.0 * Units.A, 2.0 * Units.m/Units.s)
  f_new = compute_force(x_new) # we presume compute_force() automatically attaches units to the current force
  while (t < tmax): # comparison uses natural internal units -- no need to convert
    # Update old/new
    (x_old, f_old, v_old) = (x_new, f_new, v_new)
  
    # Update positions.  All conversions are effectively automatic!
    x_new = x_old + v_old * dt + 0.5 * dt**2 * f_old / m
    # Update force - compute_force() attaches units to the computed force.
    f_new = compute_force(x_new)
    # Update velocities
    v_new = v_old + (f_old + f_new) / (2.0 * m) * dt

    # Increment integrated trajectory length.
    t += dt

SUPPORTED UNITS

Supported units and abbreviations are as follows:

  SI PREFIXES
  ato, femto, pico, nano, micro, milli, centi, deci, deca, hecto, kilo, mega, giga, tera, peta

  SIMPLE UNITS
  angle: radians (rad), degrees (deg)
  length: meters (m, cm, mm, um, nm, pm, fm), Angstroms (A), inch, foot (ft), yard (yd)
  mass: grams (kg, g), atomic mass units (u, amu, Da), stone
  time: seconds (s, ms, us, ns, ps, fs)
  charge: Coulomb (C), electronic charge units (e)
  number: mole (mol)
  temperature: kelvin (K)

  COMPOUND UNITS
  volume: liter (l, ml, ul, nl, pl)
  concentration: molar (M, mM, uM, nM, pM), molal  
  force: Newton (N), dyne, pound (lb, lbs)
  pressure: pascal (Pa, kPa, MPa, GPa), psi, atmospheres (atm), torr, bar, millimeters of mercury (mmHg)
  energy: joule (J, kJ), calorie (cal, kcal), erg
  power: watt (W)
  electric potential: volt (V), statvolt
  electric current: ampere (amp)
  magnetic field: tesla (T), gauss, oersted  
  resistance: ohm
  capacitance: farad (F)
  inductance: henry

IMPLEMENTATION DETAILS

The internal representation is SI units (related to mks).
This means that each of the following units is multiplied by 1.0 if the unit of measure is equal to:

  length: meter
  mass: kilogram
  time: second
  charge: Coulomb
  angle: radian

Many compound units come from a useful table in the back cover of Purcell "Electricity and Magnetism", Second Edition.
"""

#=============================================================================================
# Necessary fundamental constants
#=============================================================================================

# pi
from math import pi

# Avogadro's number
N_A = 6.022045e23 
Na = N_A # alias

#=============================================================================================
# SI prefixes
#=============================================================================================

ato   = 1.0e-18
femto = 1.0e-15
pico  = 1.0e-12
nano  = 1.0e-9
micro = 1.0e-6
milli = 1.0e-3
centi = 1.0e-2
deci  = 1.0e-1
deca  = 1.0e1
hecto = 1.0e2
kilo  = 1.0e3
mega  = 1.0e6
giga  = 1.0e9
tera  = 1.0e12
peta  = 1.0e15

#=============================================================================================
# SIMPLE UNITS
#=============================================================================================

#=============================================================================================
# Angles
#=============================================================================================

# Fundamental unit.
radian = 1.0 # radian
radians = radian # plural of radian
rad = radians # abbreviation of radians

# Commonly-used units.
degree = pi / 180.0 * radians # one degree, in radians
degrees = degree # plural of degree
deg = degrees # abbreviation of degrees

#=============================================================================================
# Length
#=============================================================================================

# Fundamental unit
meter = 1.0 # meter
metre = meter # English spelling of meter
m = meter # common abbreviation for meter

# Commonly-used SI unit abbreviations
cm = centi * meter # centimeter
mm = milli * meter # millimeter
um = micro * meter # micrometer
nm = nano * meter # nanometer
pm = pico * meter # picometer
fm = femto * meter # femtometer

# Other common units
Angstrom = 1.0e-10 * meter # Angstrom
Angstroms = Angstrom # plural of Angstrom
A = Angstrom # common abbreviation for Angstrom

# English units
inch = 2.54 * cm # inch
inches = inch # plural of inch
foot = 12. * inch # foot
feet = foot # plural of foot
ft = foot # common abbreviation for foot
yard = 3.0 * feet # yard
yd = yard # abbreviation for yard 

#=============================================================================================
# Mass
#=============================================================================================

# Fundmamental unit
gram = 1.0 # gram
g = gram # common abbreviation for gram

# Commonly-used SI unit abbreviations
kg = kilo * gram  # kilogram
mg = milli * gram # milligram
ng = nano * gram # nanogram
pg = pico * gram # picogram
fg = femto * gram # femtogram

# atomic mass unit
amu = (1.0 / N_A) * gram # atomic mass unit
u = amu # common abbreviation for atomic mass unit
dalton = amu # Daltons
Da = dalton # common abbreviation for dalton
kDa = kilo * dalton # kilodaltons
MDa = mega * dalton # megadaltons

# Other crazy units
stone = 6.35029318 * kg # stone in Imperial system

#=============================================================================================
# Time
#=============================================================================================

# Fundamental unit
second = 1.0 # second
seconds = second # plural of second
s = second # common abbreviation for second

# Commonly-used SI unit abbreviations
ms = milli * second # millisecond
us = micro * second # microsecond
ns = nano * second # nanosecond
ps = pico * second # picosecond
fs = femto * second # femtosecond

#=============================================================================================
# Charge
#=============================================================================================

# Fundamental unit
Coulomb = 1.0 # Coulomb
C = Coulomb # common abbreviation for Coulomb

# CGS units
esu = 1.0 / (3.0e9 * Coulomb) # esu

# Conventional units
e = 1.6021892e-19 * Coulomb # elementary (electronic) charge

#=============================================================================================
# Quantity
#=============================================================================================

mole = N_A # mole
moles = mole # plural of mole
mol = mole # common abbreviation for mole

mmol = milli * moles # millimoles
umol = micro * moles # micromoles
nmol = nano * moles # nanomoles
pmol = pico * moles # picomoles
fmol = femto * moles # femtomoles

#=============================================================================================
# Temperature
#=============================================================================================

# Fundamental unit: Kelvin
# Note that no other temperature scales are supported because they are not multiplicative rescalings.
# This unit is provided for syntactical consistency, rather than really being useful in conversions.
kelvin = 1.0 # Kelvin
Kelvin = kelvin
K = kelvin

#=============================================================================================
# COMPOUND UNITS
#=============================================================================================

#=============================================================================================
# Volume
#=============================================================================================

# SI units
liter = 1.0e3 * cm**3 # liter
litre = liter # English spelling of liter
l = liter # common abbreviation for liter
ml = milli * liter # milliliter
ul = micro * liter # microliter
pl = pico * liter # picoliter

#=============================================================================================
# Concentration
#=============================================================================================

molar = moles / liter # molar
M = molar # common abbreviation for molar

# Commonly-used SI abbreviations
mM = milli * molar # millimolar
uM = micro * molar # micromolar
nM = nano * molar # nanomolar
pM = pico * molar # picomolar
fM = femto * molar # femtomolar

molal = moles / kg # molal

#=============================================================================================
# Force
#=============================================================================================

# SI units
newton = kg * m / s**2 # newton
Newton = newton
newtons = newton
Newtons = newton
N = Newton

# Commonly uses SI abbreviations.
pN = pico * newton # piconewtons

# CGS units
dyne = g * cm / s**2 # dyne

# English units
pound = 9.80665 * m / s**2 # pound (gravity at sea level)
pounds = pound
lb = pound # common abbreviation for pound
lbs = lb # common plural abbreviation for pound

#=============================================================================================
# Pressure
#=============================================================================================

# Pascal
pascal = Newton / meter**2
Pa = pascal # pascal
kPa = kilo * pascal # kilopascal
MPa = mega * pascal # megapascal
GPa = giga * pascal # gigapascal

psi = lbs / inch**2 # pounds per square inch
atmospheres = 101325 * pascal # atmospheres
atm = atmospheres
torr = atm / 760.0 # Torr
mmHg = torr # mm of Hg
bar = 750.062 * torr
millibar = milli * bar

#=============================================================================================
# Energy
#=============================================================================================

# SI units
joule = kg * meter * second**2 # joule
Joule = joule # capitalization of joule
joules = joule # plural of joule
Joules = joule # plural of joule
J = joule # commonly-used abbreviation of joule
kJ = kilo * J # kilojoule

# CGS units
erg = dyne * cm # erg
ergs = erg

# Common units
calorie = 4.184 * joule # calorie
calories = calorie
cal = calorie # commonly used abbreviation of calorie
kcal = kilo * calorie # kilocalories
Cal = kcal # commonly used by FDA

#=============================================================================================
# Power
#=============================================================================================

watt = joule / second # watt
Watt = watt # capitalization of watt
watts = watt # plural of watt
Watts = watt
W = watt # commonly-used abbreviation for watt

#=============================================================================================
# Electric potential
#=============================================================================================

# SI units
volt = J / Coulomb # volt
Volt = volt
V = volt # commonly-used abbreviation for volt

# CGS units
statvolt = erg / esu

#=============================================================================================
# Electric current
#=============================================================================================

ampere = Coulomb / second # ampere
amperes = ampere # plural of ampere
amp = ampere # commonly-used abbreviation for ampere
A = ampere # commonly-used abbreviation for ampere
Ampere = ampere

#=============================================================================================
# Magnetic field
#=============================================================================================

# SI unit
tesla = Newton / meter / ampere # tesla
Tesla = tesla # capitalization of tesla
T = tesla # commonly-used abbreviation for tesla

# CGS unit
gauss = dyne / esu # gauss
Gauss = gauss # capitalization of gauss
oersted = dyne / esu # oersted
Oersted = oersted # capitalization of oersted

#=============================================================================================
# Resistance
#=============================================================================================

ohm = volt / amp
ohms = ohm
Omega = ohm

#=============================================================================================
# Capacitance
#=============================================================================================

farad = ohm / second
farads = farad
Farad = farad
Farads = farad
F = farad

#=============================================================================================
# Inductance
#=============================================================================================

henry = ohm * second
Henry = henry
henries = henry
Henries = henry
H = henry
