#!/usr/bin/env python

Usage = """Makes a pymol image of a conformation.

Usage     : imageconf.py [--rotate] [--cartoon] [--color=X] PDBFILE1 PDBFILE2 ...
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()

import commands, sys, os, tempfile
import scripttools
from string import *
from numpy import *
import protein, geometry

PYMOL_PATH = os.environ.get('PYMOL_PATH', '/home6/voelzv/local/bin/pymol -c')

Args = scripttools.ParseArgs(sys.argv)
PdbFiles = [fn for fn in Args["ARGS"] if not "imageconf" in fn.lower() 
            and not "python" in fn.lower()]
print "\n".join(["Could not find %s" % fn for fn in PdbFiles
                 if not os.path.isfile(fn)])
PdbFiles = [fn for fn in PdbFiles if os.path.isfile(fn)]
DoPrep = "rotate" in Args["FLAGS"]
JustCartoon = "cartoon" in Args["FLAGS"]
Color = Args.get("color", None)



def PrepPdb(Pdb):
  """Rotates principal axes to be in viewing screen."""
  p = protein.ProteinClass(Pdb=Pdb)
  #first rotate so smallest eigenvector is pointing to screen
  Pos = p.ResPos()
  evals, evecs = linalg.eig(cov(Pos, rowvar=False, bias=True))
  v = evecs[:, argmin(abs(evals))]
  Vec, Ang = geometry.GetVecMapping([0.,0.,1.], v)
  p.Rotate(Vec, Ang)
  #next rotate so larget eigenvector is pointing horozontal
  Pos = p.ResPos()
  evals, evecs = linalg.eig(cov(Pos, rowvar=False, bias=True))
  v = evecs[:, argmax(abs(evals))]
  Vec, Ang = geometry.GetVecMapping([1.,0.,0.], v)
  p.Rotate(Vec, Ang)
  #determine if we need to flip
  Pos = p.ResPos()
  avgz = average(Pos[:,2])
  frac = (Pos[:,2] > avgz).astype(float).sum() / len(Pos)
  if frac > 0.5:
    p.Rotate([1.,0.,0.], 180.)
  #save
  f, fn = tempfile.mkstemp(suffix=".pdb")
  f = os.fdopen(f, "w")
  f.close()
  p.WritePdb(fn)
  return fn



if len(PdbFiles) == 0: sys.exit()

ftmp, pymolinfile = tempfile.mkstemp(suffix=".pml")
ftmp = os.fdopen(ftmp, "w")
print "Making pictures of " + ", ".join(PdbFiles)
tempfiles = [pymolinfile]
for filename in PdbFiles:
    outfile = replace(os.path.basename(filename), '.pdb', '.png')
    if DoPrep:
      filename = PrepPdb(filename)
      tempfiles.append(filename)
    filebase = replace(os.path.basename(filename),'.pdb','')
    ftmp.write( 'load '+filename+'\n')
    if JustCartoon:
      ftmp.write( 'cmd.hide("all")\n' )
    ftmp.write( 'cmd.show(\"cartoon\"   ,\"'+filebase+'\")\n' )
    ftmp.write( 'bg_color white\n' )
    ftmp.write( 'cmd.set(\'\'\'depth_cue\'\'\',\'\'\'0\'\'\',quiet=0)\n' )
    if Color is None:
      ftmp.write( 'cmd.color("gray","all")\n' )
      ftmp.write( 'select hphobic, (resn ALA)|(resn MET)|(resn PHE)|(resn LEU)|(resn CYS)|(resn ILE)|(resn VAL)|(resn TRP)|(resn TYR)\n' )
      ftmp.write( 'cmd.color("green","hphobic")\n' )
      ftmp.write( 'select pos, (resn LYS)|(resn HIS)|(resn ARG)\n' )
      ftmp.write( 'cmd.color("blue","pos")\n' )
      ftmp.write( 'select neg, (resn GLU)|(resn ASP)\n' )
      ftmp.write( 'cmd.color("red","neg")\n' )
    elif Color == "rainbow" or Color == "spectrum":
      ftmp.write( 'cmd.spectrum()\n' )
    else:
      ftmp.write( 'cmd.color("%s","all")\n' % Color.lower() )
    ftmp.write( 'select sidechains, not (name ca,n,c,o)\n' )
    ftmp.write( 'util.cbag("sidechains")\n' )
    ftmp.write( 'deselect\n' )
    ftmp.write( 'set ray_opaque_background, off\n' )

    ftmp.write( 'ray\npng '+outfile+'\n' )
    ftmp.write( 'delete '+filebase+'\n' )
    ftmp.write( 'delete hphobic\n' )
    ftmp.write( 'delete pos\n' )
    ftmp.write( 'delete neg\n' )
    ftmp.write( 'delete sidechains\n' )
ftmp.write( 'quit\n' )
ftmp.close()

stdout = commands.getoutput(PYMOL_PATH+' '+pymolinfile)
# print stdout
for fn in tempfiles:
  os.remove(fn)

