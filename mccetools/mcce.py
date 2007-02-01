import tempfile
import os
import commands
import re
import math

def paramgen(pdbpath,mcce_location):
	"""Generate a dictionary of parameters used by MCCE.
		
		parameters:
			pdbpath - The path to the PDB file used as input for MCCE.
			mcce_location - The base directory in which MCCE is installed.
		
		returns:
			A dictionary containing a list of default MCCE parameters.
		
		"""
	# This was just generated from the sample .prm files in the
	# MCCE distribution
	mcce_params={}
	mcce_params['MCCE_HOME']=mcce_location
	mcce_params['INPDB']=pdbpath
	mcce_params['DO_PREMCCE']='t'
	mcce_params['DO_ROTAMERS']='t'
	mcce_params['DO_ENERGY']='t'
	mcce_params['DO_MONTE']='t'
	mcce_params['EXTRA']=mcce_params['MCCE_HOME']+'/extra.tpl'
	mcce_params['RENAME_RULES']=mcce_params['MCCE_HOME']+'/name.txt'
	mcce_params['EPSILON_PROT']='8.0'
	mcce_params['TITR_TYPE']='ph'
	mcce_params['TITR_PH0']='7.0'
	mcce_params['TITR_PHD']='1.0'
	mcce_params['TITR_EH0']='0.0'
	mcce_params['TITR_EHD']='30.0'
	mcce_params['TITR_STEPS']='1'
	mcce_params['MINIMIZE_SIZE']='t'
	mcce_params['TERMINALS']='t'
	mcce_params['CLASH_DISTANCE']='2.0'
	mcce_params['H2O_SASCUTOFF']='0.05'
	mcce_params['ROT_SPECIF']='f'
	mcce_params['ROT_SWAP']='t'
	mcce_params['PACK']='f'
	mcce_params['ROTATIONS']='6'
	mcce_params['SAS_CUTOFF']='1.00'
	mcce_params['VDW_CUTOFF']='10.0'
	mcce_params['REPACKS']='5000'
	mcce_params['REPACK_CUTOFF']='0.01'
	mcce_params['HDIRECTED']='f'
	mcce_params['HDIRDIFF']='1.0'
	mcce_params['HDIRLIMT']='36'
	mcce_params['HV_RELAX_NCYCLE']='0'
	mcce_params['HV_RELAX_DT']='10'
	mcce_params['HV_RELAX_NITER']='10'
	mcce_params['HV_RELAX_VDW_THR']='2.'
	mcce_params['HV_RELAX_HV_VDW_THR']='5.'
	mcce_params['HV_TORS_SCALE']='20.'
	mcce_params['HV_RELAX_N_SHAKE']='10000'
	mcce_params['HV_RELAX_CONSTRAINT']='1.0'
	mcce_params['HV_RELAX_CONSTRAINT_FRC']='20.'
	mcce_params['RELAX_H']='t'
	mcce_params['RELAX_E_THR']='-5.0'
	mcce_params['RELAX_NSTATES']='50'
	mcce_params['RELAX_N_HYD']='6'
	mcce_params['RELAX_CLASH_THR']='5.'
	mcce_params['RELAX_PHI']='1.0'
	mcce_params['RELAX_NITER']='300'
	mcce_params['RELAX_TORQ_THR']='0.5'
	mcce_params['EPSILON_SOLV']='80.0'
	mcce_params['GRIDS_DELPHI']='65'
	mcce_params['GRIDS_PER_ANG']='2.0'
	mcce_params['RADIUS_PROBE']='1.4'
	mcce_params['IONRAD']='2.0'
	mcce_params['SALT']='0.15'
	mcce_params['DELPHI_EXE']=mcce_params['MCCE_HOME']+'/bin/delphi'
	mcce_params['BIG_PAIRWISE']='5.0'
	mcce_params['MONTE_SEED']='-1'
	mcce_params['MONTE_T']='298.15'
	mcce_params['MONTE_FLIPS']='3'
	mcce_params['MONTE_NSTART']='100'
	mcce_params['MONTE_NEQ']='300'
	mcce_params['MONTE_REDUCE']='0.001'
	mcce_params['MONTE_RUNS']='6'
	mcce_params['MONTE_NITER']='2000'
	mcce_params['MONTE_TRACE']='50000'
	mcce_params['NSTATE_MAX']='1000000'
	mcce_params['DELPHI_START']='1'
	mcce_params['DELPHI_END']='99999'
	return mcce_params

def print_prm(mcceparams):
	"""Test method to dump mcceparams dictionary to screen"""
	for i in mcceparams.keys():
		print mcceparams[i]+" ("+i+")"

def run_mcce(mcceparams):
	"""Runs MCCE in the current directory, using run parameters
		given in variable mcceparams (from the paramgen function).
		This method will clobber file run.prm in the current directory,
		as well as any previous MCCE output files.
		It will generate all sorts of output files, which it is the caller's
		responsibility to clean up."""
	runfile=open("run.prm","wt")
	for i in mcceparams.keys():
		runfile.write(mcceparams[i]+" ("+i+")\n")
	runfile.close()
	output=commands.getoutput(mcceparams['MCCE_HOME']+"/bin/mcce")
	os.unlink("run.prm")
	return output

def ps_mostlikely(fort38path):
	"""Finds the set of most likely protonation states from
		the file fort38path/fort.38. Assumes one set of pH values
		in the given file.
		Returns a 4-tuple with the following elts:
			RE string for most likely states
			RE string for neutral states
			RE string for positive states
			RE string for negative states"""
	
	f=open(fort38path+"/fort.38","rt")
	lines=f.readlines()[1:]		#Skip the header
	lines2=lines[:] # Make a copy of file contents because the extraction procedure will destroy this copy.
	f.close()
	stack=[]
	# We want any _000 line - they're common to all protonation states
	searchstring="_000 "
	
	# Overview - grab a line and compare to the last one
	# if it's the same residue, keep the more likely one
	# if a new residue, add the stored line to the master "keep" list
	# and store what we just found
	while (len(lines)>0):
		nextline=lines.pop(0)
		if (len(stack)==0):
			stack.append(nextline)
		else:
			stackline=stack.pop()
			nextident=nextline[0:3]+nextline[6:10]
			stackident=stackline[0:3]+stackline[6:10]
			if (nextident!=stackident):
				# Add the identifier to the list of atoms to pull
				searchstring=searchstring+"|"+stackline[0:3]+r" "+stackline[5:15]
				#print "Keeping "+stackline[0:15]
				stack.append(nextline)
			else:
				# Next option for the same residue
				stackprob=float(stackline.split(" ")[1])
				nextprob=float(nextline.split(" ")[1])
				if (nextprob > stackprob):
					stack.append(nextline)
				else:
					stack.append(stackline)
	# Handle the last element on the stack
	stackline=stack.pop()
	searchstring=searchstring+"|"+stackline[0:3]+r" "+stackline[5:15]

	#Generate the REs for neutral/pos/neg states based on charge labeling in our duplicate array
	# Select the lines for any particular state, then transform each line into an RE, then join lines with | operator
	neutrallines = "_000 |"+("|".join(map(lambda x: x[0:3]+r" "+x[5:15],filter(lambda x:re.search(r'^...0',x),lines2))))
	poslines = "|".join(map(lambda x: x[0:3]+r" "+x[5:15],filter(lambda x:re.search(r'^...\+',x),lines2)))
	neglines = "|".join(map(lambda x: x[0:3]+r" "+x[5:15],filter(lambda x:re.search(r'^...-',x),lines2)))

	return (searchstring,neutrallines,poslines,neglines)
	
	
def renumber_atoms(pdbarr):
	"""Renumbers atom entries (list items) of a PDB file so they are
		in sequence starting from 1
		Mutates parameter pdbarr."""
	for i in range(len(pdbarr)):
		pdbarr[i]=pdbarr[i][0:6]+str(i+1).rjust(5)+pdbarr[i][11:]

def rename_termini(npdb):
	"""Renames the 'NTR' and 'CTR' residues generated by MCCE
		Renames the O and OXT atoms at the C-terminus to OC1 and OC2
		for AMBER"""
	ntrname = npdb[1][0][17:20]
	ctrname = npdb[-2][0][17:20]
	for j in [0,1]:
		for i in range(len(npdb[j])):
			npdb[j][i] = npdb[j][i][0:17]+'N'+ntrname+npdb[j][i][21:]
	for j in [-1,-2]:
		for i in range(len(npdb[j])):
			npdb[j][i] = npdb[j][i][0:17]+'C'+ctrname+npdb[j][i][21:]
		
	
def nest_pdb(pdbarr):
	nestedpdb=[]
	residue=[]
	for line in pdbarr:
		if (len(residue) == 0):
			residue.append(line)
		else:
			if (line[17:27] == residue[-1][17:27]):
				residue.append(line)
			else:
				nestedpdb.append(residue)
				residue=[line]
	nestedpdb.append(residue)
	return nestedpdb

def unnest_pdb(npdb):
	pdbarr=[]
	for res in npdb:
		for atm in res:
			pdbarr.append(atm)
	return pdbarr

def rename_residues(pdbarr):
	"""Stub to handle renaming of residues into AMBER format
		Returns new pdbarr"""
	print "At start of renaming, pdbarr has",len(pdbarr),"items"
	npdb=nest_pdb(pdbarr)
	print "Nest/unnest leads to",len(unnest_pdb(npdb)),"items"
	rename_charged(npdb)
	histidine_search(npdb)
	disulfide_search(npdb)
	rename_termini(npdb)
	pdb2=unnest_pdb(npdb)
	print "At end of renaming, pdbarr has ",len(pdb2),"items"
	return unnest_pdb(npdb)

def rename_charged(npdb):
	"""
	Generate AMBER-specific residue names for charged residues.
	
	parameters:
		npdb - "nested" PDB data structure generated by nest_pdb. Modified to reflect AMBER names.
	"""
	
	resname_and_state = map(lambda x:x[-1][17:20]+x[-1][-1],npdb)
	#print resname_and_state
	for i in range(len(npdb)):
		if (resname_and_state[i]=='HIS+'):
			npdb[i]=map(lambda x:x.replace('HIS','HIP'),npdb[i])
		elif (resname_and_state[i]=='LYS0'):
			npdb[i]=map(lambda x:x.replace('LYS','LYN'),npdb[i])
		elif (resname_and_state[i]=='LYS+'):
			npdb[i]=map(lambda x:x.replace('LYS','LYP'),npdb[i])
		elif (resname_and_state[i]=='CYS0'):
			npdb[i]=map(lambda x:x.replace('CYS','CYN'),npdb[i])
		elif (resname_and_state[i]=='CYS-'):
			npdb[i]=map(lambda x:x.replace('CYS','CYM'),npdb[i])
		elif (resname_and_state[i]=='ASP0'):
			npdb[i]=map(lambda x:x.replace('ASP','ASH'),npdb[i])
		#Aspartate requires no sub
		elif (resname_and_state[i]=='GLU0'):
			npdb[i]=map(lambda x:x.replace('GLU','GLH'),npdb[i])
		# Glutamate requires no sub
	return None

def get_coords(atomname,residue):
	for n in range(len(residue)):
		if (residue[n][12:16].strip() == atomname.strip()):
			iX=float(residue[n][30:38])
			iY=float(residue[n][38:46])
			iZ=float(residue[n][46:54])
			return (iX,iY,iZ)
	raise "Atom not found!"
	

def disulfide_search(npdb, min_dist = 1.8, max_dist = 2.2):
	"""rename CYS to CYX in the presence of disulfide bonds
	1.8-2.2 ang"""
	residues_to_rename=set([])
	for i in range(len(npdb)):
		if (npdb[i][0][17:20] != 'CYS'):
			continue
		# Found a cysteine, now track down the sulfur
		iX,iY,iZ=get_coords('SG',npdb[i])
		
		for j in range(i+1,len(npdb)):
			if (npdb[j][0][17:20]!= 'CYS'):
				continue
			jX,jY,jZ=get_coords('SG',npdb[j])

			dX=iX-jX
			dY=iY-jY
			dZ=iZ-jZ
			distance = math.sqrt(dX*dX+dY*dY+dZ*dZ)
			if (distance >= min_dist and distance <= max_dist):
				residues_to_rename = residues_to_rename | set([i,j])
	
	# Rename the residues we selected
	for i in residues_to_rename:
		npdb[i]=map(lambda x:x.replace('CYS','CYS2'),npdb[i])		
	return None

def histidine_search(npdb):
	"""Tries to rename HIS to HID or HIE
		Assumes that there are hydrogens in this structure"""
	for i in range(len(npdb)):
		resname = npdb[i][0][17:20]
		if (resname == 'HIS'):
			epsilonh_found=0
			for atom in npdb[i]:
				if (atom[13:16] == 'HE2'):
					epsilonh_found = 1
			if (epsilonh_found):
				npdb[i]=map(lambda x:x.replace('HIS','HIE'),npdb[i])
			else:
				npdb[i]=map(lambda x:x.replace('HIS','HID'),npdb[i])
	return None

def pdb_cleanup(pdbarr):
	"""handle removal of extraneous entries from
		MCCE PDB files
		Returns new pdbarr"""
	# Nuke everything after the coordinates
	# Use a default occupancy of 1 and a default B-factor of 0
	for i in range(len(pdbarr)):
		atom=pdbarr[i][0:54]
		# The atom symbol will be the first character after stripping
		# whitespace and numbers
		# NOTE: this will fail for atoms with more than one letter in their symbol
		#            eg, counterions
		atomsymbol = atom[12:16].strip(" 0123456789")[0]
		atom = atom + "1.00".rjust(6) + "0.00".rjust(6) + " "*10 + atomsymbol.rjust(2)+ " "*2
		pdbarr[i]=atom
	npdb=nest_pdb(pdbarr)
	for i in range(len(npdb)):
		# Set the residue index numbers
		# NOTE: We set chain identifier to blank, so this won't work for multi-chain PDBs
		for j in range(len(npdb[i])):
			atom = npdb[i][j]
			atom = atom[:21]+ " " + str(i+1).rjust(4)+ " "*4 + atom[30:]
			#print len(atom)
			npdb[i][j] = atom
		
	return unnest_pdb(npdb)
	
def protonation_state(pdbfile,ph,mccepath,cleanup=True):
	"""Finds the ML protonation state of a given PDB file at a
		given pH by running MCCE.
		Returns: list containing lines of the output PDB file"""
	# Set up the parameters for the MCCE run
	params=paramgen(pdbfile,mccepath)
	params['TITR_PH0']=str(ph)
	params['TITR_STEPS']="1"
	
	# Create a temporary directory and run MCCE
	cwd=os.getcwd()
	tempdir=tempfile.mkdtemp();
	print tempdir
	os.chdir(tempdir)
	run_mcce(params)
	
	pdbarr = ps_processmcce(tempdir)

	# Generate a breakpoint for interactive testing...
#	raw_input("about to clean house...")
	
	# clean up the temp dir
	# There should only be files, no subdirs
	if (cleanup):
		for i in os.listdir(tempdir):
			os.unlink(i)
		os.chdir(cwd)
		os.rmdir(tempdir)
	
	return pdbarr

def ps_processmcce(tempdir):
	"""Handles the file processing work for protonation_state"""
	# Build and use a regex to grab the appropriate entries from the MCCE PDB
	sstr,neut,pos,neg=ps_mostlikely(tempdir)
	rx=re.compile(sstr)
	print "searchstring: \"",sstr
	f=open(tempdir+"/step2_out.pdb","rt")
	pdbarr=filter(lambda x:rx.search(x),f.readlines())
	f.close()
	
	# Label the lines of the PDB file with charge state
	rxn=re.compile(neut)
	rxp=re.compile(pos)
	rxm=re.compile(neg)
	for i in range(len(pdbarr)):
		if (rxn.search(pdbarr[i])):
			pdbarr[i]=pdbarr[i]+"0"
		elif (rxp.search(pdbarr[i])):
			pdbarr[i]=pdbarr[i]+"+"
		elif (rxm.search(pdbarr[i])):
			pdbarr[i]=pdbarr[i]+"-"
		else:
			print "ERROR!"
			raise "Missed a case in charge state labeling"
	
	renumber_atoms(pdbarr)
	pdbarr=rename_residues(pdbarr)
	pdbarr=pdb_cleanup(pdbarr)
	return pdbarr
