#!usr/bin/python

import sys


def read_grofile( grofile ):
    
    atomnums = {}
    resnames = {}
    fin = open(grofile,'r')
    lines = fin.readlines()
    lines.pop(0)
    lines.pop(0)
    lines.pop()
    for line in lines:
	resnum = int(line[0:5])
	resname = line[5:9].strip()
        fields = line.split()
        if len(fields) > 1:
	    # print fields
	    atomname = fields[1]
	    atomnum = int(fields[2])
	    atomnums[ (resnum,atomname) ] = atomnum
	    resnames[ resnum ] = resname
    fin.close()
    return atomnums, resnames


def read_noefile( noefile ):
    
    noepairs = []
    fin = open(noefile,'r')
    lines = fin.readlines()
    for line in lines:
      if line[0] != '#':
        fields = line.split()
        if len(fields) > 1:
	    resname1   = fields[0]
	    resnum1    = int(fields[1])
	    atomname1  = fields[2]
	    resname2   = fields[3]
	    resnum2    = int(fields[4])
	    atomname2  = fields[5]
	    strength = fields[6]
	    noepairs.append( [(resnum1, atomname1), (resnum2, atomname2), strength] )
    fin.close()
    return noepairs

def PairsToRest( atomnums, noepairs, resnames):
    
    gro_atompairs = []	# atomnum1, atomnum2, strength
    pymol_atompairs = []    # resi1, name1, resi2, name2
    xml_atompairs = []    # resi1, name1, resi2, name2
    for pair in noepairs:
	found1 = False
	found2 = False
	pair1 = pair[0]
	pair2 = pair[1]
	oldpair1 = None
	oldpair2 = None
	if atomnums.has_key( pair1 ):
            found1 = True
	else:
            oldpair1 = pair1
	    pair1 = (pair[0][0], pair[0][1]+'1')
            if atomnums.has_key( pair1 ):
		found1 = True
	
	if atomnums.has_key( pair2 ):
            found2 = True
	else:
            oldpair2 = pair2
            pair2 = (pair[1][0], pair[1][1]+'1')
            if atomnums.has_key( pair2 ):
        	    found2 = True
	
	if (found1 & found2):
            atom1 = atomnums[ pair1 ]
            atom2 = atomnums[ pair2 ]
            strength = pair[2]
            gro_atompairs.append( (atom1, atom2, strength) )
	    pymol_atompairs.append( (pair1[0], pair1[1], pair2[0], pair2[1]) )   # resi1, name1, resi2, name2
	    tmp = (pair1[0], resnames[pair1[0]], pair1[1], pair2[0], resnames[pair2[0]], pair2[1])
	    xml_atompairs.append( tmp )   # resi1, name1, resi2, name2
            
            # print '--->    HIT!', (atom1, atom2, strength)
            # print '    pair1',oldpair1, '->', pair1
            # print '    pair2', oldpair2, '->', pair2
	else:
            print '-> can\'t find:'
            print '     pair1',oldpair1, '->', pair1
            print '     pair2', oldpair2, '->', pair2
    
    print 'gro_atompairs', gro_atompairs
    print 'len(gro_atompairs)', len(gro_atompairs)
    
    print 'pymol_atompairs', pymol_atompairs
    print 'len(pymol_atompairs)', len(pymol_atompairs)
    
    print 'xml_atompairs', xml_atompairs
    print 'len(xml_atompairs)', len(xml_atompairs)
    
    return [gro_atompairs, pymol_atompairs, xml_atompairs]


def WriteRestFile(gro_atompairs, outfile):
    
    outstr =          '[ distance_restraints ]\n'
    outstr = outstr + '; ai   aj    type  index  type\'  low   up1  up2   fac\n'
    index = 0
    for p in gro_atompairs:
       index=index+1
       up1 = 0.18
       up2 = 1.0
       if p[2] == 'Very':
           up2 = 0.6
       if p[2] == 'Weak':
           up2 = 0.5
       if p[2] == 'Medium':
           up2 = 0.33
       if p[2] == 'Strong':
           up2 = 0.27
       outstr=outstr +'%-6d %-6d %-6d %-6d %-6d %2.2f  %4.2f %4.2f %4.2f\n'%(p[0],p[1],1,index,1,0.0,up1,up2,1.0)
    
    fout = open(outfile,'w')
    fout.write(outstr)
    fout.close()

def WritePyMOLFile(pymol_atompairs, outfile):
    
    fout = open(outfile,'w')
    index = 0
    for p in pymol_atompairs:
       index=index+1
       fout.write( 'select resi%d_%s, resi %d and name %s\n'%(p[0],p[1],p[0],p[1]) )
       fout.write( 'select resi%d_%s, resi %d and name %s\n'%(p[2],p[3],p[2],p[3]) )
       fout.write( 'distance d%d = resi%d_%s, resi%d_%s\n'%(index,p[0],p[1],p[2],p[3]) )
    fout.close()

def WriteXMLFile(xml_atompairs, outfile):
		    
    outstr = '[ distance_restraints ]\n'
    outstr = outstr + '; ai   aj    type  index  type\'  low   up1  up2   fac\n'
    index = 0
    for p in gro_atompairs:
       index=index+1
       up1 = 0.18
       up2 = 1.0
       if p[2] == 'Very':
           up2 = 0.6
       if p[2] == 'Weak':
           up2 = 0.5
       if p[2] == 'Medium':
           up2 = 0.33
       if p[2] == 'Strong':
           up2 = 0.27
       outstr=outstr +'%-6d %-6d %-6d %-6d %-6d %2.2f  %4.2f %4.2f %4.2f\n'%(p[0],p[1],1,index,1,0.0,up1,up2,1.0)

    fout = open(outfile,'w')
    fout.write(outstr)
    fout.close()



if __name__ == '__main__':
  
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("-g", "--grofile", dest="grofile", help="the *.gro input file", metavar="FILE")
  parser.add_option("-n", "--noefile", dest="noefile", help="a text file with the NOEs listed in the following format: \"NARG 1  HG      THR 12  HB      Weak\nPHE 2   HA\"      These names must correspond to the *.gro input file.""", metavar="FILE")
  parser.add_option("-o", "--outname", dest="outname", help="the *.regress input file", metavar="FILE")
  parser.add_option("-v", "--verbose", dest="verbose", help="verbose", action="store_true", default=False )
  (options, args) = parser.parse_args()
  
  restfile = options.outname + '.itp'
  pymolfile = options.outname + '.log'
  xmlfile = options.outname + '.xml'
  
  atomnums, resnames = read_grofile(options.grofile)
  print 'atomnums',atomnums
  
  keys = atomnums.keys()
  keys.sort()
  for key in keys:
      print key, atomnums[key]
  
  noepairs = read_noefile(options.noefile)
  print 'len(noepairs)',len(noepairs)
  
  # convert the pairs to GROMACS .top restraints
  [gro_atompairs, pymol_atompairs, xml_atompairs] = PairsToRest( atomnums, noepairs, resnames)
  
  WriteRestFile(gro_atompairs, restfile)
  WritePyMOLFile(pymol_atompairs, pymolfile)
  WriteXMLFile(xml_atompairs, xmlfile)




