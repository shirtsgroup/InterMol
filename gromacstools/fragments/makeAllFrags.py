#!/usr/bin/env python

import sys, os, glob, commands, string
from zam import protein

# needed for tlc2olc
import mmtools.pdbtools as pdbtools
from mmtools.gromacstools.TopologyFile import *

usage = """makeAllFrags.py OUTDIR [TYPE]

    Will build PDB (*.pdb) and gmx topology files (*.itp) for every amino acid fragment

    OUTDIR    The output directory

    TYPE      "sidechain" - Replaces all backbone atoms as extra C_beta hydrogen.  Default. 
           or "backbone"  - caps N- and C- termini with neutral hydrogen.
"""

# constants

olc2tlc = {'A':'ALA',
           'C':'CYS',
           'D':'ASP',
           'E':'GLU',
           'F':'PHE',
           'G':'GLY',
           'H':'HIS',
           'I':'ILE',
           'K':'LYS',
           'L':'LEU',
           'M':'MET',
           'N':'ASN',
           'P':'PRO',
           'R':'ARG',
           'S':'SER',
           'T':'THR',
           'V':'VAL',
           'W':'TRP',
           'Y':'TYR'}



# functions


def MakeFragment(seq, Outdir, Name, FragType='sidechain', forcefield='amber03', Verbose=False, CleanUp=True):
    """Makes a fragment PDB and topoolgy files for a given sequence,
    with asdfghj
    
    REQUIRED
    
    Seq         one-letter amino acid seqiuence
    Outdir      Directory to write output files
    Name        basename for filenames
    
    PARAMETERS
    
    FragType    'sidechain' or 'backbone'
    Verbose     Show extended print statements
    Cleanup     If True, delete pesky "#*" gromacs backup files
    """
    
    # make outdir if it doesn't exist yet
    if not os.path.exists(Outdir):
        os.mkdir(Outdir)
        
    p = protein.ProteinClass(Seq=seq)
    pdbfile = os.path.join(Outdir, Name+'.pdb')
    p.WritePdb(pdbfile)
    
    
    
    # Process the PDB file to make in more palatable to pdb2gmx
    
    # 1) delete all hydrogrens except backbone 'H'
    #    NOTE: the zam hydrogen names don't necessarily match the gmx names, so let pdb2gmx fill them in correctly
    # 2) change LYS to LYP
    # 3) change CYS to CYN
    # 4) change HIS to HIP
    PdbAtoms = pdbtools.readAtomsFromPDB(pdbfile)
    i = 0
    while i < len(PdbAtoms):
        if PdbAtoms[i]["resName"] == 'LYS ':
            PdbAtoms[i]["resName"] = 'LYP '
        if PdbAtoms[i]["resName"] == 'CYS ':
            PdbAtoms[i]["resName"] = 'CYN '
        if PdbAtoms[i]["resName"] == 'HIS ':
            PdbAtoms[i]["resName"] = 'HIP '
        AtomName = PdbAtoms[i]["name"].strip()
        if (AtomName.count('H') > 0) & (len(AtomName) > 1) & (AtomName.count('NH') == 0) \
            & (AtomName.count('CH') == 0) & (AtomName.count('OH') == 0) :   # leave backbone H
            PdbAtoms.pop(i)
        else:
            i+=1
        pdbtools.writeAtomsToPDB(pdbfile, PdbAtoms, renumber=True)
    
    
    # make a gro and top file from the pdbfile
    grofile = pdbfile.replace('.pdb','.gro')
    topfile = pdbfile.replace('.pdb','.top')
    cmd = 'pdb2gmx -ff %s -f %s -o %s -p %s'%(forcefield, pdbfile, grofile, topfile)
    print '>>',cmd
    output = commands.getoutput(cmd)
    print output
    # clean up backup files
    backupfiles = glob.glob(os.path.join(Outdir,'#*'))
    for bfile in backupfiles:
        os.remove(bfile)

    # read in a topology file object
    mytop = TopologyFile(topfile)
    
    # read in a topology file object
    mygro = pdbtools.GromacsStructureFromGrofile(grofile)
    
    
    AtomLookup = {}
    for atomnum in mytop.molecules[0].atoms.keys():
        AtomLookup[ mytop.molecules[0].get_atomname(atomnum) ] = atomnum
    
        
        
        
    if FragType=='sidechain':

        BackboneAtoms = ['N','H','HA','C','O']
        print 'BackboneAtoms',BackboneAtoms
        
        # Remove the backbone atoms from the grofile

        for AtomName in BackboneAtoms:
            remove_atomnum = None
            if Verbose: print 'Looking for AtomName =', AtomName
            for atom in mygro.atoms:
                if atom.atomname.strip() == AtomName:
                    if Verbose: print 'Found one!  atom',atom.atomnum,'is', AtomName
                    remove_atomnum = atom.atomnum
                    
                    if Verbose: print 'Removing', AtomName, 'atomnum =', remove_atomnum
                    mygro.removeatom(remove_atomnum-1) #ACK! <-- this is in list numbering! but the grofile numbering starts at 1.
            mygro.write(grofile)
            mygro = pdbtools.GromacsStructureFromGrofile(grofile)
            if Verbose: print mygro

        if Verbose: print '%%% Final grofile after deleting backbone %%%'
        if Verbose: print mygro
        
        
        # change 'CA' to a 'HB' atom in both topfile and grofile
        
        ### Lookup hydrogen parameters from atom 'HB1'
        if AtomLookup.has_key('HB1'):
            newHydrogenAtomType = mytop.molecules[0].get_atomtype(AtomLookup['HB1'])
            newHydrogenAtomMass = mytop.molecules[0].get_mass(AtomLookup['HB1'])
        elif AtomLookup.has_key('HB'):
            newHydrogenAtomType = mytop.molecules[0].get_atomtype(AtomLookup['HB'])
            newHydrogenAtomMass = mytop.molecules[0].get_mass(AtomLookup['HB'])
            # change the 'HB' name to 'HB1', becase we're going to add one more.
            mytop.molecules[0].set_atomname(AtomLookup['HB'],'HB1')
            for i in range(len(mygro.atoms)):
                if mygro.atoms[i].atomname.strip() == 'HB':
                    mygro.atoms[i].atomname = '%6s'%'HB1'
                    mygro.write(grofile)
                    mygro = pdbtools.GromacsStructureFromGrofile(grofile)
                    break

            
        else:
            raise "Can't find any 'HB1' or 'HB' atom in the topology file!"
        
        ### find out if there's a free 'HB#' name
        count = 2
        newHydrogenName = 'HB%d'%count
        while AtomLookup.has_key(newHydrogenName):
            count += 1
            newHydrogenName = 'HB%d'%count

        ### change the 'CA' atom to the newHydrogen in the topfile
        if Verbose: print 'Changing CA atom to', newHydrogenName
        AtomLookup[newHydrogenName] = AtomLookup['CA'] # add the new hydrogen to the lookup
        if not AtomLookup.has_key('CA'):
            raise "Can't find atom 'CA' in the topology file!"
        mytop.molecules[0].set_atomname(AtomLookup['CA'], newHydrogenName)
        mytop.molecules[0].set_atomtype(AtomLookup['CA'], newHydrogenAtomType)
        mytop.molecules[0].set_mass(AtomLookup['CA'], newHydrogenAtomMass)
        # don't change the charge *yet*, as we need to rebalance the charge over all the 'HB#' atoms
        
        ### change 'CA' to 'HB#' in the grofile
        if Verbose: print "Changing 'CA' to 'HB#' in the grofile..."
        for i in range(len(mygro.atoms)):
            if mygro.atoms[i].atomname.strip() == 'CA':
                if Verbose: print 'Found one!  atom',mygro.atoms[i].atomnum,'is', mygro.atoms[i].atomname
                mygro.atoms[i].atomname = '%6s'%newHydrogenName
                mygro.write(grofile)
                mygro = pdbtools.GromacsStructureFromGrofile(grofile)
                break
            
        print "%%% Final grofile after Changing 'CA' to 'HB#' in the grofile... %%%"
        print mygro
        
        # compile the total charge on all the 'HB#' atoms and the backbone
        HBAtoms = []
        for AtomName in AtomLookup.keys():
            if AtomName.count('HB')>0: HBAtoms.append(AtomName)    
        if Verbose: print 'HBAtoms',HBAtoms
        BackboneCharges = [ mytop.molecules[0].get_charge(AtomLookup[AtomName]) for AtomName in BackboneAtoms ]
        if Verbose: print 'BackboneCharges',BackboneCharges
        HBCharges = [ mytop.molecules[0].get_charge(AtomLookup[AtomName]) for AtomName in HBAtoms ]
        if Verbose: print 'HBCharges',HBCharges
        ChargePerHBAtom = (sum(BackboneCharges) + sum(HBCharges))/len(HBAtoms)
        if Verbose: print 'ChargePerHBAtom', ChargePerHBAtom

        # set all the HB partial charges to the new ChargePerHBAtom value
        for AtomName in HBAtoms:
            mytop.molecules[0].set_charge(AtomLookup[AtomName], ChargePerHBAtom)
        
        # delete the backbone atoms
        for AtomName in BackboneAtoms:
            mytop.molecules[0].remove_atom(AtomLookup[AtomName])

        # Renumber the atoms according the to grofile
        mytop.molecules[0].renumber_atoms(NameOrdering=[(1,atom.atomname.strip()) for atom in mygro.atoms])
        
        print '%%% NEW TOPOLOGY FILE %%%'
        print mytop
        
        mygro.write(grofile)
        mytop.write(topfile)
    
    
    
    
    if FragType=='backbone':
        
        # build a tpr for the grofile (so we can protonate it next)
        
        minimization_mdpfile = """title               =  minimization
cpp                 =  cpp
integrator          =  steep
nsteps              =  10000

; Energy minimizing stuff
emtol               =  100.0
emstep              =  0.001

comm_mode           =  none
nstcomm             =  1
nstxout             =  500
nstvout             =  0
nstfout             =  0
nstlog              =  50000
nstenergy           =  50000
nstxtcout           =  50000
xtc_grps            =  Protein
nstlist             =  10
ns_type             =  grid
pbc                 =  xyz  ; means yes
coulombtype         =  PME
epsilon_r           =  80
vdwtype             =  shift
; Berendsen temperature coupling is on in two groups
Tcoupl              =  berendsen
tc-grps             =  Protein
tau_t               =  0.1         
ref_t               =  308   
; Energy monitoring
energygrps          =  Protein
; Isotropic pressure coupling is now on
;Pcoupl             =  no
Pcoupltype          =  isotropic
tau_p               =  1
compressibility     =  4.5e-5
ref_p               =  1.01325

; Generate velocites is on at 300 K.
gen_vel             =  no
gen_temp            =  308
gen_seed            =  173529

; Duan-Kollman used hbonds not all-bonds
; but we wanted to increase dt (earlier)
; so I've always used all-bonds

constraints = none
"""
        
        fout = open('minimize.mdp','w')
        fout.write(minimization_mdpfile)
        fout.close()
        
        # get rid of 'spc.itp' #include lines 
        fin = open(topfile,'r')
        toplines = fin.readlines()
        fin.close()
        fout = open(topfile,'w')
        for line in toplines:
            if line.count('spc') == 0:
                fout.write(line)
        fout.close()
        
        # edit the grofile so it has a bigger box
        cmd = 'editconf -f %s -box 5.0 5.0 5.0 -o %s'%(grofile,grofile)
        print '>>', cmd
        output = commands.getoutput(cmd)
        print output
        
        # build a tprfile
        tprfile = grofile.replace('.gro','.tpr')
        cmd = 'grompp -f minimize.mdp -c %s -p %s -o %s'%(grofile,topfile,tprfile)
        print '>>', cmd
        output = commands.getoutput(cmd)
        print output
        
        # protonate the tprfile
        grofile_protonated = grofile.replace('.gro','_protonated.gro')
        cmd = 'echo "0\n" | protonate -s %s -o %s'%(tprfile,grofile_protonated)
        print '>>', cmd
        output = commands.getoutput(cmd)
        print output
        

        # NOTE: The "protonate" command gives a *.gro file ONLY (no topfile)
        #       It adds a H to the amino N, and adds an O- to the carboyl carbon
        # We need to:
        #       1) add the new H's to the topology file
        #       2) change the new O to a H in the grofile; get rid of the N-termial H3
        
        example = """Example protonated grofile
   22
    1LEU      N    1   2.554   2.227   2.418
    1LEU     H1    2   2.459   2.195   2.418
    1LEU     H2    3   2.601   2.193   2.336
    1LEU     H3    4   2.601   2.193   2.500
    1LEU     CA    5   2.556   2.372   2.418
    1LEU     HA    6   2.508   2.397   2.334
    1LEU     CB    7   2.486   2.428   2.541
    1LEU    HB1    8   2.532   2.396   2.624
    1LEU    HB2    9   2.391   2.397   2.542
    1LEU     CG   10   2.490   2.580   2.537
    1LEU     HG   11   2.585   2.611   2.537
    1LEU    CD1   12   2.419   2.630   2.412
    1LEU   HD11   13   2.422   2.730   2.409
    1LEU   HD12   14   2.465   2.594   2.331
    1LEU   HD13   15   2.324   2.599   2.412
    1LEU    CD2   16   2.419   2.636   2.661
    1LEU   HD21   17   2.422   2.736   2.658
    1LEU   HD22   18   2.324   2.605   2.661
    1LEU   HD23   19   2.465   2.604   2.743
    1LEU      C   20   2.699   2.425   2.418
    1LEU     O1   21   2.776   2.361   2.418
    1LEU     O2   22   2.716   2.524   2.418
   5.00000   5.00000   5.00000"""

        # sys.exit(1)
        
        mygro = pdbtools.GromacsStructureFromGrofile(grofile_protonated)
        
        # get rid of 'H3' on the N-terminal residue
        for i in range(len(mygro.atoms)):
            if mygro.atoms[i].resname.strip() != 'PRO':
                if (mygro.atoms[i].atomname.strip() == 'H3') & (mygro.atoms[i].resnum == 1):
                    mygro.removeatom(i)
                    break
            else:
                if (mygro.atoms[i].atomname.strip() == 'H2') & (mygro.atoms[i].resnum == 1):
                    mygro.removeatom(i)
                    break
        mygro.write(grofile)

        print '%%% AFTER H3 removal %%%'
        print mygro

        
        # Change 'O1' on the C-terminal residue to O
        # Change 'O2' on the C-terminal residue to H
        
        ### find last resnum
        resnums = []
        for i in range(len(mygro.atoms)):
            if resnums.count(mygro.atoms[i].resnum) == 0:
                resnums.append(mygro.atoms[i].resnum)
        resnums.sort()
        lastres = resnums[-1]
        
        for i in range(len(mygro.atoms)):
            if (mygro.atoms[i].atomname.strip() == 'O1') & (mygro.atoms[i].resnum == lastres):
                mygro.atoms[i].atomname = '   O'
            if (mygro.atoms[i].atomname.strip() == 'O2') & (mygro.atoms[i].resnum == lastres):
                mygro.atoms[i].atomname = '  H4'
        mygro.write(grofile)
        
        print '%%% AFTER O1 -> O and O2 -> H  %%%'
        print mygro
        
        
        
        
        ResNums = []
        for atomnum in mytop.molecules[0].atoms.keys():
            ThisResNum = mytop.molecules[0].get_atomresnum(atomnum)
            if ResNums.count(ThisResNum)==0:
                ResNums.append(ThisResNum)
        ResNums.sort()
        FirstResNum = ResNums[0]
        LastResNum = ResNums[-1]
        
        
        # Add the N-terminal hydrogen to the topology file
        if mygro.atoms[0].resname.strip() == 'PRO':
            newHydrogenName = 'H1'
            newHydrogenAtomType = 'amber99_17'
            newHydrogenAtomMass = 1.008
        else:
            newHydrogenName = 'H2'
            for atomnum in mytop.molecules[0].atoms.keys():
              if mytop.molecules[0].get_atomresnum(atomnum) == FirstResNum:
                if mytop.molecules[0].get_atomname(atomnum) == 'H1':
                    newHydrogenAtomType = mytop.molecules[0].get_atomtype(atomnum)
                    newHydrogenAtomMass = mytop.molecules[0].get_mass(atomnum)
                if mytop.molecules[0].get_atomname(atomnum) == 'H':
                    newHydrogenAtomType = mytop.molecules[0].get_atomtype(atomnum)
                    newHydrogenAtomMass = mytop.molecules[0].get_mass(atomnum)
                    mytop.molecules[0].set_atomname(atomnum,'H1')

        newAtomNum = mytop.molecules[0].get_new_atomnum()
        mytop.molecules[0].atoms[newAtomNum] = copy.copy(mytop.molecules[0].atoms[newAtomNum-1])
        mytop.molecules[0].set_atomname(newAtomNum, newHydrogenName)
        mytop.molecules[0].set_atomtype(newAtomNum, newHydrogenAtomType)
        mytop.molecules[0].set_cgnr(newAtomNum, newAtomNum)
        mytop.molecules[0].set_atomnum(newAtomNum, newAtomNum)
        mytop.molecules[0].set_atomresnum(newAtomNum, FirstResNum)
        mytop.molecules[0].set_mass(newAtomNum, newHydrogenAtomMass)
        mytop.molecules[0].set_charge(newAtomNum, 0.0)
        
        ### make a lookup for the neighbor atoms we need to connect
        Neighbor = {}   # dict to be filled as {'N': atomnum}, e.g.
        for atomnum, atomlist in mytop.molecules[0].atoms.iteritems():
            if mytop.molecules[0].get_atomresnum(atomnum) == FirstResNum:
                if mytop.molecules[0].get_atomname(atomnum) == 'N':
                    Neighbor['N'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'C':
                    Neighbor['C'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'CA':
                    Neighbor['CA'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'CB':
                    Neighbor['CB'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'HA':
                    Neighbor['HA'] = atomnum
        
        ### connect bonds
        mytop.molecules[0].bonds[( Neighbor['N'], newAtomNum) ] = [ str(Neighbor['N']), str(newAtomNum), '1']
        ### connect pairs
        mytop.molecules[0].pairs[( newAtomNum, Neighbor['CB']) ] = [ str(newAtomNum), str(Neighbor['CB']), '1']
        mytop.molecules[0].pairs[( newAtomNum, Neighbor['HA']) ] = [ str(newAtomNum), str(Neighbor['HA']), '1']
        mytop.molecules[0].pairs[( newAtomNum, Neighbor['C']) ] = [ str(newAtomNum), str(Neighbor['C']), '1']
        ### connect angles
        print 'Adding angle:', ( newAtomNum, Neighbor['N'], Neighbor['CA'])
        mytop.molecules[0].angles[( newAtomNum, Neighbor['N'], Neighbor['CA']) ] = [ str(newAtomNum), str(Neighbor['N']), str(Neighbor['CA']), '1']
        ### connect dihedrals
        mytop.molecules[0].dihedrals[( newAtomNum, Neighbor['N'], Neighbor['CA'], Neighbor['CB']) ] = [ str(newAtomNum), str(Neighbor['N']), str(Neighbor['CA']), str(Neighbor['CB']), '3']
        mytop.molecules[0].dihedrals[( newAtomNum, Neighbor['N'], Neighbor['CA'], Neighbor['HA']) ] = [ str(newAtomNum), str(Neighbor['N']), str(Neighbor['CA']), str(Neighbor['HA']), '3']
        mytop.molecules[0].dihedrals[( newAtomNum, Neighbor['N'], Neighbor['CA'], Neighbor['C']) ] = [ str(newAtomNum), str(Neighbor['N']), str(Neighbor['CA']), str(Neighbor['C']), '3']
        
        
            
        print '%%% AFTER first (N-term) H added to topfile %%%'
        print mytop
            
        # Add the C-terminal hydrogen to the topology file
        newHydrogenName = 'H4'
        newHydrogenAtomType = 'amber99_23'
        newHydrogenAtomMass = 1.008
        
        newAtomNum = mytop.molecules[0].get_new_atomnum()
        mytop.molecules[0].atoms[newAtomNum] = copy.copy(mytop.molecules[0].atoms[newAtomNum-1])
        mytop.molecules[0].set_atomname(newAtomNum, newHydrogenName)
        mytop.molecules[0].set_atomtype(newAtomNum, newHydrogenAtomType)
        mytop.molecules[0].set_cgnr(newAtomNum, newAtomNum)
        mytop.molecules[0].set_atomnum(newAtomNum, newAtomNum)
        mytop.molecules[0].set_atomresnum(newAtomNum, LastResNum)
        mytop.molecules[0].set_mass(newAtomNum, newHydrogenAtomMass)
        mytop.molecules[0].set_charge(newAtomNum, 0.0)
        
        ### make a lookup for the neighbor atoms we need to connect
        Neighbor = {}   # dict to be filled as {'C': atomnum}, e.g.
        for atomnum, atomlist in mytop.molecules[0].atoms.iteritems():
            if mytop.molecules[0].get_atomresnum(atomnum) == LastResNum:
                if mytop.molecules[0].get_atomname(atomnum) == 'N':
                    Neighbor['N'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'C':
                    Neighbor['C'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'CA':
                    Neighbor['CA'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'CB':
                    Neighbor['CB'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'HA':
                    Neighbor['HA'] = atomnum
                if mytop.molecules[0].get_atomname(atomnum) == 'O':
                    Neighbor['O'] = atomnum
        
        ### connect bonds
        mytop.molecules[0].bonds[( Neighbor['C'], newAtomNum) ] = [ str(Neighbor['C']), str(newAtomNum), '1']
        ### connect pairs
        mytop.molecules[0].pairs[( newAtomNum, Neighbor['N']) ] = [ str(newAtomNum), str(Neighbor['N']), '1']
        mytop.molecules[0].pairs[( newAtomNum, Neighbor['CB']) ] = [ str(newAtomNum), str(Neighbor['CB']), '1']
        mytop.molecules[0].pairs[( newAtomNum, Neighbor['HA']) ] = [ str(newAtomNum), str(Neighbor['HA']), '1']
        ### connect angles
        print 'Adding angle:', ( newAtomNum, Neighbor['C'], Neighbor['O'])
        mytop.molecules[0].angles[( newAtomNum, Neighbor['C'], Neighbor['O']) ] = [ str(newAtomNum), str(Neighbor['C']), str(Neighbor['O']), '1']
        print 'Adding angle:', ( newAtomNum, Neighbor['C'], Neighbor['CA'])
        mytop.molecules[0].angles[( newAtomNum, Neighbor['C'], Neighbor['CA']) ] = [ str(newAtomNum), str(Neighbor['C']), str(Neighbor['CA']), '1']
        ### connect dihedrals
        mytop.molecules[0].dihedrals[( newAtomNum, Neighbor['C'], Neighbor['CA'], Neighbor['N']) ] = [ str(newAtomNum), str(Neighbor['C']), str(Neighbor['CA']), str(Neighbor['N']), '3']
        mytop.molecules[0].dihedrals[( newAtomNum, Neighbor['C'], Neighbor['CA'], Neighbor['CB']) ] = [ str(newAtomNum), str(Neighbor['C']), str(Neighbor['CA']), str(Neighbor['CB']), '3']
        mytop.molecules[0].dihedrals[( newAtomNum, Neighbor['C'], Neighbor['CA'], Neighbor['HA']) ] = [ str(newAtomNum), str(Neighbor['C']), str(Neighbor['CA']), str(Neighbor['HA']), '3']

        
        # Renumber the atoms according the to grofile
        mytop.molecules[0].renumber_atoms(NameOrdering=[(atom.resnum, atom.atomname.strip()) for atom in mygro.atoms])
 
                
        print '%%% AFTER second (C-term) H added to topfile %%%'
        print mytop
        
        mygro.write(grofile)
        mytop.write(topfile)
        
        
    ################
    # MINIMIZATION #
    ################

    
    # get rid of 'spc.itp' #include lines 
    fin = open(topfile,'r')
    toplines = fin.readlines()
    fin.close()
    fout = open(topfile,'w')
    for line in toplines:
        if line.count('spc') == 0:
            fout.write(line)
    fout.close()

    # edit the grofile so it has a bigger box
    cmd = 'editconf -f %s -box 5.0 5.0 5.0 -o %s'%(grofile,grofile)
    print '>>', cmd
    output = commands.getoutput(cmd)
    print output
    
    # build a tprfile for minimization
    tprfile = grofile.replace('.gro','.tpr')
    cmd = 'grompp -f minimize.mdp -c %s -p %s -o %s'%(grofile,topfile,tprfile)
    print '>>', cmd
    output = commands.getoutput(cmd)
    print output

    # minimize the structure
    grofile_minimized = grofile.replace('.gro','_min.gro')
    cmd = 'mdrun -s %s -c %s'%(tprfile,grofile_minimized)
    print '>>', cmd
    output = commands.getoutput(cmd)
    print output

    
    
        
        
                
    # Cleanup
    backupfiles = glob.glob('./#*') + glob.glob( os.path.join(Outdir,'#*') )
    if Verbose: print '### CLEANUP ###'
    for bfile in backupfiles:
        if Verbose: print 'removing', bfile
        os.remove(bfile)
        
        
        

    
    

if __name__ == '__main__':
    
    
    if len(sys.argv) < 2:
        print usage
        sys.exit()
        
    outdir = sys.argv[1]
    fragtype = 'sidechain'
    if len(sys.argv) >= 3:
        fragtype = sys.argv[2]
         
    Verbose = True
    
    if Verbose:
        print 'outdir', outdir
        print 'fragtype', fragtype
        
    sequences = copy.copy(olc2tlc.keys())
    if (fragtype == 'sidechain'):
        sequences.remove('P')
        sequences.remove('G') # proline and glycine won't work with fragments
        
    for seq in ['W']:   # for testing
    # for seq in sequences:  
        
        if (fragtype == 'sidechain'):
            if len(seq) > 1:
                print "Only single-amino acid sequences are allowed with fragtype = %s.  Skipping..."%fragtype
                continue
            if seq.count('P') > 0:
                print "Proline is not allowed with fragtype = %s.  Skipping..."%fragtype
                continue
            if seq.count('G') > 0:
                print "GLycine is not allowed with fragtype = %s.  Skipping..."%fragtype
                continue
        MakeFragment(seq, outdir, seq, FragType=fragtype)
        
        