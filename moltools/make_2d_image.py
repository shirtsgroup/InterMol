#!/usr/bin/python

#Version 1 (rough) D. Mobley 6/25/2007

from openeye.oechem import *
from openeye.oedepict import *
from openeye.oeiupac import *
import os

def name_to_img(name, imgfile):
  """Take a IUPAC name and write out an image file containing a 2d structure of that molecule, currently in GIF format, to the specified imgfile output."""

  #Generate molecule
  molecule = OEMol()
  status = OEParseIUPACName(molecule, name)

  #This loosely based on mol2gif.py in the wrappers/python/examples directory

  #Set up various stuff (not too sure what it does)
  OESetDimensionFromCoords(molecule)
  OESuppressHydrogens(molecule)
  OEDepictCoordinates(molecule)
  OEMDLPerceiveBondStereo(molecule)

  #initialize view
  view = OEDepictView()
  #Adjust size of bonds/scale
  view.SetPixScale(56)
  #Turn off display of logo
  view.SetLogo(False)
  view.SetAromaticDashes(True)
  view.SetMolecule(molecule)
  img = OE8BitImage()
  view.RenderImage(img, True, 0, 0)

  ofs = oeosstream()
  #For some reason I can't currently get PDF output to come out without giving a big blank background, some I am using gifs for now.
  #OEWriteEPS(ofs, img)
  OEWriteGIF(ofs, img, -1)
  s = ofs.str()
  #open(imgfile,'w').write(s)
  open(imgfile,'w').write(s)
  #os.system('ps2pdf test.eps')

