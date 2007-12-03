#=============================================================================================
# gencoil.py 
#
# Generates a random-coil PDB structure from an amino acid sequence .
#=============================================================================================
# AUTHORS
#
# Written by Vincent Voelz, Stanford University, 2007-07-06
#=============================================================================================
# REQUIREMENTS
#=============================================================================================
# TODO
#=============================================================================================
# VERSION CONTROL INFORMATION
__version__ = "$Revision: $"                                                                         
                     
#=============================================================================================
# IMPORTS

import sys, os
import random, math, string, tempfile 
import numpy

from optparse import OptionParser    # For parsing of command line arguments

from mmtools.utilities import Units
from mmtools.pipelinetools.RamaTable import *
from mmtools.pipelinetools.DihedralEnergy import *

from mmtools.pdbtools import *
from mmtools.gromacstools.GromacsData import *


# We need ZAM functions to modify dihedral angles, etc. 
sys.path.append('../zam')
import protein, proteinfunc, proteinconst, pdbtools, coords

 
#=============================================================================================
# FUNCTIONS
#=============================================================================================


def maximumDistance( PDBfile ):
    """Return the maximum interatomic distance of a PDB file."""
    pdbcoords = coords.GetPdbCoords( PDBfile )
    distances = getAllPairwiseDistances( pdbcoords )
    if distances != None:
        return math.sqrt( distances.max() )
    else:
        return None 


def getAllPairwiseDistances( coords ):
    """Use the law of cosines to get the distances bweteen all atoms:    
    d_{ij}^2 = x_i^2 + y_j^2 - 2\dot(x_i, x_j)    
    Returns a matrix D of squared-distances D_ij.
    """
    
    origin = coords.sum(axis=0)/float(coords.shape[0])       
    # print 'origin', origin
    coords = coords - origin
    # The incoming coords, X, are an N x 3 matrix.  
    transpose_coords = numpy.transpose(coords)    # X^T
    outer = numpy.dot(coords, transpose_coords)   # X(X^T)
                
    squaredDistances = numpy.zeros( outer.shape )

    dotprodrows = numpy.dot( coords**2, numpy.ones( transpose_coords.shape ) )   # M_ij = x_i*x_i
    dotprodcols = numpy.transpose(dotprodrows)      # M_ij = x_j*x_j
    squaredDistances = dotprodrows + dotprodcols - outer
    # print 'squaredDistances', squaredDistances
    return squaredDistances


def gencoil(sequence, outdir, nconf=10, temperature=300.0, dimension=None, useNN='single', useLib='native'):
    """Generate a series of [nconfs] unfolded conformations in directory outdir.

    ARGUMENTS
    sequence            one-letter sequence string, either a filename, or a string
    outdir		The name of the directory where outfiles files are written.
   
    OPTIONS
    nconfs		The number (int) of unfolded conformations to generate.
    temperature         The temperature of the ensemble (T_ref=300K for dihedral statistical potential)
    dimension           The (float) maximum atom-to-atom distance allowed for an unfolded conformation, in Angstroms.
                        This is useful to generating unfolded conformations in a finite box
    useNN               If 'double' (default), account for nearest-neighbor (NN) effects when generating dihedrals,
                        following Jha, Sosnick, et al. PNAS (2005).   If useNN='single', ignore these correlations.
                        
    useLib              'coil'    - Dihedral random coil library 
                        'C_A'
                        'C_AB'
                        'C_ABT'
                        'native'
                       
                         
    """

    if dimension != None:
        dimension = float(dimension)

    kB = 1.987 * Units.calorie    # (...per mol per degree)  
    kT = kB * temperature * Units.K / Units.kcal
    kT_ref = kB * 300. * Units.K / Units.kcal    # our reference is 300K
    print 'effective kT =', kT,     '(%3.1f)'%temperature
    print 'reference kT =', kT_ref, '(300.0 K)'

    if os.path.exists(sequence):
        fseq = open(sequence,'r')
        seqstr = fseq.readline().strip()
        fseq.close()
    else:
        seqstr = sequence

    glycineString = ''
    for i in range(0,len(seqstr)):
        glycineString = glycineString+'G'

    # Set MC parameters
    Verbose = False
    PrintEvery = 1         # print progress to screen every PrintEvery iterations
    WritePDBEvery = 10     # Write a PDB file every WritePDBEvery snapshots 
    EquilIters = 50      # start writing PDB files after EquilIters
    AngleRange = 30.0  # Range of angles to take phi-psi MC moves from
    NTweak = 3         # Change NTweak residues' dihedral angles at a time 
    WritePDBEvery = len(seqstr)/NTweak     # Write a PDB file every WritePDBEvery snapshots 

    
    print 'Running %d iterations to equilibrate before accepting viable conformations...'%EquilIters
    print 'Conformations of atoms less than 0.65 Angstroms are not accepted.'
    print 'Angle range to take phi-psi MC moves from: %f'%AngleRange

    # Prepare output directory files
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    tmpPDBfile = os.path.join(outdir,'tmp.pdb')
    fene = open(os.path.join(outdir,'out.ene'),'w')
    eneheader = '%-16s %-16s %-16s %-16s %-8s\n'%('iteration','Energy','num_overlaps','dimension','wrotePDB')
    fene.write(eneheader)

    # Create a ZAM Protein Object
    p = protein.ProteinClass(Seq=glycineString)

    # Create a DihedralEnergy Object 
    ene = DihedralEnergy(useLib=useLib)
  
    # Create a RamaTable objects
    rama = RamaTable()

    # mess up the angles to start with -- draw from uncorrelated (no NN) distribution)
    for i in range(0,len(p)):
        Phi, Psi = p.PhiPsi(i)
        if not Phi is None and not Psi is None:
            Phi, Psi, prob = rama.getRandomDihedral()
            print 'initial prob for residue',i,'=',prob
            p.RotateToPhiPsi(i, Phi, Psi)

    npdb = 0
    iters = 0
    while npdb < nconf:
        
        
        # Each iteration, correct the center of mass
        centerPos = (p.Pos.sum(0)/float(p.Pos.shape[0]))
        p.Pos = p.Pos - centerPos
        
        
        iters = iters + 1
        oldEnergy = 0.
        newEnergy = 0.

        if useNN=='single':
          # Propose a move of all the dihedral angles
          oldDihedrals = []
          newDihedrals = []
          # remember the old dihedrals, and calculate the old energy
          for i in range(0,len(p)):
            Phi, Psi = p.PhiPsi(i)
            rescode = string.capitalize(seqstr[i])
            oldDihedrals.append( (Phi,Psi) )
            if not Phi is None and not Psi is None:
                basin = rama.basinFromDihedral(Phi, Psi)
                oldEnergy = oldEnergy + ene.single(rescode, basin) 
          # propose a change in the dihedrals
          # METHOD ONE:  Pick N random residues at a time, and tweak them.
          for i in range(0, len(oldDihedrals)):
              newDihedrals.append( oldDihedrals[i] )
          for n in range(0, NTweak):
              ipick = random.randint(0,len(p)-1)
              Phi, Psi, prob = rama.getRandomDihedral( basin=rama.getRandomBasin() )
              newDihedrals[ipick] = (Phi,Psi) 
          # calculate the energy of the new dihedrals
          for i in range(0,len(newDihedrals)):
              (Phi, Psi) = newDihedrals[i]
              basin = rama.basinFromDihedral(Phi, Psi)
              newEnergy = newEnergy + ene.single(rescode, basin) 

        elif useNN=='double':
          # Propose a move of all the dihedral angles
          oldDihedrals = []
          newDihedrals = []
          # remember the old dihedrals
          for i in range(0,len(p)-1):
            Phi1, Psi1 = p.PhiPsi(i)
            Phi2, Psi2 = p.PhiPsi(i+1)
            oldDihedrals.append( (Phi1,Psi1) )
          oldDihedrals.append( (Phi2,Psi2) )
          # remember the energy of the old dihedrals
          for i in range(0,len(oldDihedrals)-1):
            (Phi1, Psi1) = oldDihedrals[i]
            (Phi2, Psi2) = oldDihedrals[i+1]
            rescode1 = seqstr[i]
            rescode2 = seqstr[i+1]
            if not Phi1 is None and not Psi1 is None and not Phi2 is None and not Psi2 is None:
                basin1 = rama.basinFromDihedral(Phi1, Psi1) 
                basin2 = rama.basinFromDihedral(Phi2, Psi2)
                oldEnergy = oldEnergy + ene.double( rescode1, rescode2, basin1, basin2 )
          # propose a change of the dihedrals
          # METHOD ONE:  Pick one random residue and tweak it.
          for i in range(0, len(oldDihedrals)):
              newDihedrals.append( oldDihedrals[i] )
          for n in range(0, NTweak):
              ipick = random.randint(0,len(p)-1)
              Phi, Psi, prob = rama.getRandomDihedral( basin=rama.getRandomBasin() )
              newDihedrals[ipick] = (Phi,Psi)
          # calculate the energy of the new dihedrals
          for i in range(0,len(newDihedrals)-1):
            (Phi1, Psi1) = newDihedrals[i]
            (Phi2, Psi2) = newDihedrals[i+1]
            rescode1 = seqstr[i]
            rescode2 = seqstr[i+1]
            if not Phi1 is None and not Psi1 is None and not Phi2 is None and not Psi2 is None:
                basin1 = rama.basinFromDihedral(Phi1, Psi1)
                basin2 = rama.basinFromDihedral(Phi2, Psi2)
                # print 'ene.double( %s, %s, %d, %d )'%( rescode1, rescode2, basin1, basin2 ), ene.double( rescode1, rescode2, basin1, basin2 )
                newEnergy = newEnergy + ene.double( rescode1, rescode2, basin1, basin2 ) 

        else:
           raise "useNN must be either 'single' or 'double'"
           
 
        accept = False
        acceptprob = 1.0
        if newEnergy < oldEnergy:
            accept = True
        else:
            acceptprob = math.exp((1./kT)*(oldEnergy - newEnergy))
            if random.random() < acceptprob:
                accept = True

        nclashes = None
        thislength = None
        lengthOK = False
        wrotePDB = False
        
        if ((iters-1) % 100) == 0:
            progressPDBfile = os.path.join(outdir,'progress'+string.zfill(iters,8)+'.pdb')
            print '\t\t*** writing', progressPDBfile
            p.WritePdb(progressPDBfile)

        if accept == True:
            for i in range(0,len(newDihedrals)):
                p.RotateToPhiPsi(i, newDihedrals[i][0], newDihedrals[i][1])
            p.WritePdb( tmpPDBfile )

            # a "clash"/overlap is when any two atoms are within 0.6 Angstroms (according to Scott's scripts)
            # We only will write PDB files for configurations with NO clashes, because these are minimizable!
            clashes = pdbtools.GetOverlaps( tmpPDBfile )
            nclashes = len(clashes)
            thislength = maximumDistance( tmpPDBfile )
            if Verbose: print '\t*** accepting.   nclashes =',nclashes
            if (dimension==None) or (thislength <= dimension):
                lengthOK = True
            if (nclashes == 0) and (iters >= EquilIters) and lengthOK and (iters%WritePDBEvery==0):
                outPDBfile = os.path.join(outdir,'out'+string.zfill(npdb,8)+'.pdb')
                print '\t\t*** writing', outPDBfile
                p.WritePdb(outPDBfile)
                npdb = npdb + 1
                wrotePDB = True
            fene.write( '%-16d %-16e %-16d %-16.2f %-8d\n'%(iters,newEnergy,nclashes,thislength,int(wrotePDB)) )

        # Print status of the MC search
        if iters % PrintEvery == 0:
            if accept == True:
                print 'iter',iters,':',npdb,' of',nconf,'PDB files. acceptprob','%-6e'%acceptprob,'nclashes =%-2d length = %4.2f'%(nclashes,thislength), ' | oldEnergy','%-6f'%oldEnergy,'newEnergy','%-6f'%newEnergy
            else:
                print 'iter',iters,':',npdb,' of',nconf,'PDB files. acceptprob','%-6e'%acceptprob,'NOT ACCEPTED                 | oldEnergy','%-6f'%oldEnergy,'newEnergy','%-6f'%newEnergy

    # cleanup
    fene.close()
    os.remove(tmpPDBfile)


  

    
    
    
    
    
    
    
#===========================================================================
# MAIN
#===========================================================================

if __name__ == '__main__':

    # Create command-line argument options.
    usage_string = """%prog --seq SEQUENCE --outdir OUTDIR --nconf N --temp TEMPERATURE --dimension LENGTH --useNN KEYWORD

    Generate a series of random-coil PDB conformations of a given sequence, but with
    no side chains (i.e. polyglycine). These PDB files can then be threaded using MODELLER
    e.g., to generate unfolded states with the correct sidechains.
    
    REQUIRED ARGUMENTS
    
    -s SEQUENCE              Either a protein sequence 'MKVIFLK', or a filename containing that sequence
    -o OUTDIR                The name of the output directorry.  Will create if it doesn't exist. 
    -n N                     The number of output confomrations to be written to the outdir
    -t TEMPERATURE           The temperature (K) used to the statistical scoring functio (which has reference temp 300K).
    
    OPTIONAL ARGUMENTS
    
    -d LENGTH                If specified, restrict the output conformations to those whose maximum atom-to-atom distance is
                             less than LENGTH.
    -u KEYWORD               'single' (default) for no nearest-neighbor (NN) correlations, 'double' to include them.  
"""

    version_string = "%prog %__version__"

    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-s", "--seq", metavar='SEQUENCE',
            action="store", type="string", dest='sequence', default=None,
            help="The amino acid sequence of the chain to be modeled.  \
            This can be either a one-letter code string, or a file containing that string.")
    parser.add_option("-o", "--outdir", metavar='OUTDIR',
            action="store", type="string", dest='outdir', default=None,
            help="Output directory to write the PDB conformations.  An energy *.ene file is also written.")
    parser.add_option("-n", "--nconf", metavar='N',
            action="store", type="string", dest='nconf', default=10,
            help="The number of random-coil PDB conformations to output.")
    parser.add_option("-t", "--temp", metavar='TEMPERATURE',
            action="store", type="string", dest='temperature', default=300.,
            help="The effective temperature for the MC algorithm (in Kelvin).")
    parser.add_option("-d", "--dimension", metavar='LENGTH',
            action="store", type="string", dest='dimension', default=None,
            help="Restrict conformations to those with a maximum dimension, in Angstroms. Optional.")
    parser.add_option("-u", "--useNN", metavar='KEYWORD',
            action="store", type="string", dest='useNN', default='single',
            help="'single' (default, no NN) or 'double', which considers nearest-neighbor (NN) effects.")


    # Parse command-line arguments.
    (options,args) = parser.parse_args()
    parser.get_usage()

    # Perform minimal error checking.
    if ((not options.sequence) or (not options.outdir)):
        parser.error("Must specify the sequence and the outdir.\n")
        parser.print_usage()

    # Send arguments to gencoil()
    sequence = options.sequence
    outdir = options.outdir
    nconf = int(options.nconf)
    temperature = float(options.temperature)

    gencoil( options.sequence, options.outdir, nconf=int(options.nconf), 
             temperature=float(options.temperature), dimension=options.dimension, useNN=options.useNN  )
 
