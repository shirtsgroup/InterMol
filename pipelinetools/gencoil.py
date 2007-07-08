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
# - Need a replacement PDB + dihedral manipulation scheme to replace the ZAM modules
# ONCE WE DO THIS, we can put this functionality in the pipeline module!
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
 
sys.path.append('/Users/vincentvoelz/scott/scripts')
    # For now, we need ZAM functions to modify dihedral angles, etc.
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


def gencoil(sequence, outdir, nconf=10, temperature=300.0, dimension=None):
    """Generate a series of [nconfs] unfolded conformations in directory outdir."""

    if dimension != None:
        dimension = float(dimension)

    kB = 1.987 * Units.calorie    # (...per mol per degree)  
    kT = kB * temperature * Units.K / Units.kcal
    kT_ref = kB * 300. * Units.K / Units.kcal    # our reference is 300K
    print 'effective kT =', kT
    print 'reference kT =', kT_ref

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
    PrintEvery = 50     # print progress to screen every PrintEvery iterations
    EquilIters = 100    # start writing PDB files after EquilIters

    # Prepare output directory files
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    tmpPDBfile = os.path.join(outdir,'tmp.pdb')
    fene = open(os.path.join(outdir,'out.ene'),'w')
    eneheader = '%-16s %-16s %-16s %-16s %-8s\n'%('iteration','LogProb','num_overlaps','dimension','wrotePDB')
    fene.write(eneheader)

    # Create a Protein Object
    p = protein.ProteinClass(Seq=glycineString)

    # mess up the angles to start with
    for i in range(0,len(p)):
        Phi, Psi = p.PhiPsi(i)
        if not Phi is None and not Psi is None:
            Phi += 180. * (2.*random.random() - 1.)
            Psi += 180. * (2.*random.random() - 1.)
            p.RotateToPhiPsi(i, Phi, Psi)

    npdb = 0
    iters = 0
    while npdb < nconf:

        iters = iters + 1
        oldLogProbConf = 0.
        newLogProbConf = 0.

        # Propose a move of all the dihedral angles
        oldDihedrals = []
        newDihedrals = []
        for i in range(0,len(p)):
            Phi, Psi = p.PhiPsi(i)
            oldDihedrals.append( (Phi,Psi) )
            if not Phi is None and not Psi is None:
                oldLogProbConf = oldLogProbConf + math.log( proteinfunc.GetRamaProb(Phi, Psi, proteinconst.RamaProb) )
                Phi += 180. * (2.*random.random() - 1.)
                Psi += 180. * (2.*random.random() - 1.)
                newDihedrals.append( (Phi,Psi) )
                newLogProbConf = newLogProbConf + math.log( proteinfunc.GetRamaProb(Phi, Psi, proteinconst.RamaProb) )

        if Verbose: print 'oldLogProbConf',oldLogProbConf,'newLogProbConf',newLogProbConf

        accept = False
        acceptprob = 1.0
        if newLogProbConf > oldLogProbConf:
            accept = True
        else:
            acceptprob = math.exp((kT/kT_ref)*(newLogProbConf - oldLogProbConf))
            if random.random() < acceptprob:
                accept = True

        nclashes = None
        thislength = None
        lengthOK = False
        wrotePDB = False
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
            if (nclashes == 0) and (iters >= EquilIters) and lengthOK:
                outPDBfile = os.path.join(outdir,'out'+string.zfill(npdb,8)+'.pdb')
                print '\t\t*** writing', outPDBfile
                p.WritePdb(outPDBfile)
                npdb = npdb + 1
                wrotePDB = True
            fene.write( '%-16d %-16e %-16d %-16.2f %-8d\n'%(iters,newLogProbConf,nclashes,thislength,int(wrotePDB)) )

        # Print status of the MC search
        if iters % PrintEvery == 0:
            if accept == True:
                print 'iter',iters,':',npdb,' of',nconf,'PDB files. acceptprob',acceptprob,'nclashes',nclashes,'thislength',thislength
            else:
                print 'iter',iters,':',npdb,' of',nconf,'PDB files. acceptprob',acceptprob,'NOT ACCEPTED.'

    # cleanup
    fene.close()
    os.remove(tmpPDBfile)


#===========================================================================
# MAIN
#===========================================================================

if __name__ == '__main__':

    # Create command-line argument options.
    usage_string = """%prog --seq SEQUENCE --outdir OUTDIR --nconf N --temp TEMPERATURE --dimension LENGTH

Generate a series of random-coil PDB conformations of polyglycine. These PDB
files can then be threaded using MODELLER, e.g., to generate unfolded states.

NOTES

A Monte Carlo algorithm is used  with all-residue random dihedral moves, and an effective
energy function U(phi,psi) = -kT log p(phi,psi). Each residue is considered to be independent 
(i.e. no nearest-neighbor effects).  Only conformations that have no pair of
atoms closer than 0.65 Angstroms are written to the output directory.
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


    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Perform minimal error checking.
    if ((not options.sequence) or (not options.outdir)):
        parser.error("Must specify the sequence and the outdir.\n")

    # Send arguments to gencoil()
    sequence = options.sequence
    outdir = options.outdir
    nconf = int(options.nconf)
    temperature = float(options.temperature)

    gencoil( options.sequence, options.outdir, nconf=int(options.nconf), 
             temperature=float(options.temperature), dimension=options.dimension )
 
