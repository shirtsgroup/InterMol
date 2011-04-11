from Topology.Molecule import Molecule
from Topology.System import System
from Topology.Atom import Atom
from Topology.Force import *
from Topology.Decorators import *

system = System("Test System")

m = Molecule("Molecule 1")
m2 = Molecule("Molecule 2")
m3 = Molecule("Molecule 1")

a = Atom(1, 
        "His", 
        "C",
        "12",
        1 * units.nanometers,
        2 * units.nanometers,
        3 * units.nanometers,
        1 * units.nanometers * units.picoseconds**(-1),
        2 * units.nanometers * units.picoseconds**(-1),
        3 * units.nanometers * units.picoseconds**(-1),
        1 * units.kilojoules_per_mole * units.nanometers**(-1),
        2 * units.kilojoules_per_mole * units.nanometers**(-1),
        3 * units.kilojoules_per_mole * units.nanometers**(-1))

m.addAtom(a)

force = Bond("stuff", "stuff2", 15 * units.nanometers, 10 * units.kilojoules_per_mole * units.nanometers**(-2))
print force.length

system.addMolecule(m)
system.addMolecule(m2)
system.addMolecule(m3)
print system.molecules
set =  system.molecules["Molecule 1"]
for molecule in set:
    print molecule.atoms
