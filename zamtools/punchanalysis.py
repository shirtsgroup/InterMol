#!/usr/bin/env python

#LAST MODIFIED: 05-04-07 by MSS

description = """
#
# punchanalysis.py:
#
#     Calculates the "degree fraction" i.e. the "punch" of
#     every contact in a ZAP simulation, given cooperativity
#     and stability thresholds.
#
#     NOTE:  There is NO STANDALONE/cmdline version of punchanalysis.py *yet*

"""

usage = """
Usage:  punchanalysis.py DataPath OutputPath PairList ContactCut CritCutProb CritCutCoop NFrameRead TargetTemp

    PairList            a list of tuples like [(40,52),(53,55),(43,55)] containing the residue numbers
                        of contacts to be monitored. If PairList = []
    ContactCut          the cutoff distance for the contact probabilities
    CritCutProb         probability threshold for punch calculation
    CritCutCoop         cooperativity threshold for punch calculation
    NFrameRead          the number of frames in the distance file.
    TargetTemp          the target temperature for reweighting results
"""


  
import sys
import os
import commands
import pickle
import glob
import gzip
import copy

from string  import *
import math
from random import *

from numpy import *

#sys.path.append('/home1/voelzv/md/coopanal/Config')

DEBUG = False
VERBOSE = False

GZIPPED_FILES = True       
ZAPFILES = True
MAX_SNAPSHOTS = 2000000     # this is a safegaurd to stop if the job runs too long...
PRINT_EVERY = 100
CO_LIMIT = 3                # don't consider contacts with CO < CO_LIMIT

K_BOLTZ = 1.98717e-3        # in kcal/(K*mol)
INVLOG2 = 1./math.log(2)


def RunAnal(DataPath, OutputPath, PairList, ContactCut, CritCutProb,  CritCutCoop,
                      NFrameRead, TargetTemp, ReplicaInd, LogwKN):
                      
    
    # DataPath:         This is a string like "/r41-56/data" that gives the folder where the distance files are.
    #                           Distance files are named "0.dist.txt.gz" etc.  
    # OutputPath:       This is a string like "/r41-56/anal" where you should write your results files to.
    # PairList:         This is a list of tuples like [(40,52),(53,55),(43,55)] that contains the residue
    #                           numbers of contacts to be monitored.  Actually, you probably don't need this
    #                           variable since you can just get everything from the distance file headers, but
    #                           it's helpful to have.  The number of contacts is then len(PairList).  Note that
    #                           the residue numbers here correspond to the ones in the native, full sequence AND
    #                           they start at zero (so for protein G, they run between 0 and 55).
    # ContactCut:       This is a float specifying the cutoff distance for the contact probabilities.
    # CritCutProb:      This is a float specifying the probability above which a contact is considered
    #                           punchworthy (e.g., you have been using 0.5... I thought it might be nice to have
    #                           it in the call so we can easily change this through the ZAP scripts.)
    # NFrameRead:       This integer gives the number of frames in the distance file and in LogwKN.
    # TargetTemp:       This is a float for the target temperature for reweighting results.
    # ReplicaInd:       This an array of replica indices to use.
    # LogwKN:           This is a NumPy array of floats of dimension (len(ReplicaInd),NFrameRead).  LogwKN[k,n] gives
    #                           the log of the weight of configuration n of temperature k, to a constant (which
    #                           can be removed after normalization).  This comes from WHAM.
  


    # OTHER Notes for Vince
    #
    # 1) we don't need an AtomGroup because the distances (of some pre-determined group)
    #    have been calculated ahead of time!
    # 2) files should always be GZIPPED
    # 3) Use True and False -- a la Scott

    ### check for too few pairs
    if len(PairList) < 2:
      return 0., 0., [0.]*len(PairList), [0.]*len(PairList)
    
    ### Convert the PairList (starts at 0) to residues (which start at 1)
    PairList = [(a+1,b+1) for (a,b) in PairList]

    ### create output directory if it doesn't exist ###
    if not(os.access(OutputPath,0)):                   # 0 is the dummy chmod mode...
        os.mkdir(OutputPath)

    distfilelist = glob.glob(os.path.join(DataPath, '*.dist.txt.gz'))
    if (DEBUG):
        print distfilelist

    ### get rid of any previously existing contact state dictionaries    
    dictionaryfile = os.path.join(OutputPath, 'cstates.dictionary')
    if (DEBUG):
        print 'dictionaryfile',dictionaryfile         
    if os.access(dictionaryfile,0):            # 0 is the dummy chmod mode...
        os.remove(dictionaryfile)

    contact_dictionary = {}

    ### Convert the distance trajectories into contact state trajectories
    for dtrajfile in distfilelist:

        new_dictionary = extract_states(        dtrajfile,
                                                contact_dictionary,
                                                DataPath,
                                                OutputPath,
                                                PairList,
                                                ContactCut,
                                                CritCutProb,
                                                NFrameRead,
                                                ReplicaInd )
                                                
        contact_dictionary = new_dictionary


    ########################################
    # Save the dictionary of contact states
    # ---we dont really need to do this..
    
    # fdict = open(dictionaryfile,'w')
    # pickle.dump(contact_dictionary,fdict)
    # fdict.close()
    
    ###################################################
    ### Write a contact state index to file
    ### store in 'cstates.index.txt.gz'

    indexfile = OutputPath+'/cstates.index.txt.gz'
    findex = gzip.open( indexfile, 'w')
    findex.write('index       contact-state\n')
    print 'Writing the contact state index to file', indexfile, '...'
    flipdict = {}
    for i,j in contact_dictionary.iteritems():
        flipdict[j] = i
    sorted_keys = flipdict.keys()
    sorted_keys.sort()
    for i in sorted_keys:    
        tmpstr = str(i)
        while len(tmpstr) < 12:
            tmpstr = tmpstr + ' '
        tmpstr = tmpstr + flipdict[i] + '\n'
        findex.write(tmpstr)
    findex.close() 


    ################################################
    ### calculate the free energy of each contact state and write to file
    ### store in 'cstates.free.txt.gz'

    ctrajfilelist = []
    for i in ReplicaInd:
        ctrajfilelist.append(OutputPath + '/%d.ctraj.gz' % i)
        
    maxindex = max(contact_dictionary.values())
    minindex = min(contact_dictionary.values())
    if DEBUG:
        print 'maxindex',maxindex,'minindex',minindex   
    
    Zlist = []
    for i in range(minindex, maxindex+1):
        Zlist.append( 0.0 )
    Z = array( Zlist )
    
    minLogwKN = LogwKN.min()
    weights = exp( LogwKN - minLogwKN )
    
    states_list = []
    for (k, ctrajfile) in enumerate(ctrajfilelist):
    
        states_list.append( [] )
        fctraj = gzip.open(ctrajfile,'r')
        line = fctraj.readline()
        while line != '':
            states_list[k].append( int(line) )
            line = fctraj.readline()
        fctraj.close()

    states = array( states_list )
    if DEBUG:
        print 'shape(states)',shape(states)
        print 'shape(weights)',shape(weights)

    Zindex = 0
    for state_index in range(minindex,maxindex+1):
        pickout = where(states == state_index,1,0)
        Z[Zindex] = sum(sum( weights*pickout))
        Zindex = Zindex + 1
    
    free_energies = -K_BOLTZ*TargetTemp*log(Z)
    free_energies = free_energies - free_energies.min()     
    
    freefile = OutputPath + '/cstates.free.txt.gz'
    fout = gzip.open(freefile,'w')
    fout.write('index       free-energy\n')
    for i in range(0,len(Z)):
        tmpstr = str(i)
        while len(tmpstr) < 12:
            tmpstr = tmpstr + ' '
        tmpstr = tmpstr + str(free_energies[i]) + '\n'
        fout.write(tmpstr)          
    fout.close()             

        
    ################################################
    ### calculate the cooperativity of all contact pairs
    ### store in 'coop.txt.gz'


    coop2data = pairs( free_energies, TargetTemp, freefile, indexfile, OutputPath, PairList )
    
    if DEBUG:
        print 'coop2data',coop2data
    
    PunchStabTot, PunchCoopTot, PunchStabVals, PunchCoopVals = \
      punch( coop2data, OutputPath, PairList, CritCutProb,  CritCutCoop )
    
    return [PunchStabTot, PunchCoopTot, PunchStabVals, PunchCoopVals]



def extract_states(dtrajfile, contact_dictionary, DataPath, OutputPath, PairList, ContactCut, CritCutProb,
                   NFrameRead, ReplicaInd):
                 
                 
    print 'Analyzing distance trajectory',dtrajfile,'...'

    ##########
    # Files

    if is_gzipped(dtrajfile):
        ftraj = gzip.open(dtrajfile,'r')
    else:
        ftraj = open(dtrajfile,'r')

    ### Parse out the name (i.e. number) of the trajectory
    trajstrings = splitfields(dtrajfile,'/')
    while len(trajstrings) > 1:
        trajstrings.pop(0)             ### parse out  '#.trj'
    trajstrings = splitfields(trajstrings[0],'.')

    outname = OutputPath + '/' + trajstrings[0]  ### parse out just the '#'

    outtrajfile = outname + '.ctraj.gz'

    fctraj = gzip.open(outtrajfile,'w') 


    #############################
    # Parse the trajectory

    if (DEBUG):
        print 'Parsing the distance trajectory...'

    snapshot = 0
    stopreading = 0

    header = ftraj.readline()   # the first line is junk (either \n or 'ACE' ???)
    header_fields = splitfields(header)

    contacts_by_column = []
    header_fields.pop(0)       # get rid of 'index'
    for field in header_fields:
        field = '(' + field + ')'
        contacts_by_column.append( eval( field ) )

    if (DEBUG):
        print contacts_by_column
        print PairList
    
    if len(PairList) == 0:
        PairList = contacts_by_column

    # find the max and min residue, for creating matrix distances[]
    pair =  PairList[0]
    maxres = max( pair[0], pair[1] )
    minres = min( pair[0], pair[1] )
    for pair in PairList:
        maxres = max(maxres, max(pair[0], pair[1]) )
        minres = min(minres, min(pair[0], pair[1]) )
        

    # print 'maxres',maxres,'minres',minres  
    nres = maxres - minres + 1
    distances = (ContactCut + 1.0)*ones([nres,nres], 'Float32' )
    # print 'distances',distances
    # print 'ContactCut',ContactCut

    maxkeyindex = 0         # the highest keyindex so far


    snapshot = 0
    line = ftraj.readline()
    while line != '':
        fields = splitfields(line)
        index = int(fields.pop(0))
        for i in range(0,len(fields)):
            fields[i] = float(fields[i])
            resi = contacts_by_column[i][0]
            resj = contacts_by_column[i][1]
            # print 'resi,resj',resi,resj, float(fields[i])
            distances[resi-minres,resj-minres] = float(fields[i])
            # print 'distances'  , distances
        cmatrix = distances < ContactCut
        
        # print 'cmatrix',cmatrix
        # print (cmatrix == True)
        
        
        contactstate_string = unfurl(cmatrix, minres, nres)

        tmpstr = str(getkeyindex(contactstate_string,contact_dictionary)) + '\n'

        fctraj.write(tmpstr)
        snapshot = snapshot + 1

        if snapshot % PRINT_EVERY == 0 and VERBOSE:
            print 'snapshot',snapshot,contactstate_string  
  
        line = ftraj.readline()




    ################
    # Clean up

    fctraj.close()   

    if VERBOSE: print '...Done.'
    
    return contact_dictionary




def unfurl(cmatrix, minres, nres):

    raveled_indices =  flatnonzero(cmatrix)
    cstate = []
    for i in raveled_indices:
        resi = i/nres 
        resj = i - nres*resi        
        cstate.append(   (resi+ minres,resj+ minres) )

    tmp = '['
    for entry in cstate:
        tmp = tmp + ' ' + str(entry[0]) + ' ' + str(entry[1])
    tmp = tmp + ' ]'
        
    return tmp
    

def getkeyindex(key,contact_dictionary):

    ### print 'contact_dictionary.has_key(key)',contact_dictionary.has_key(key)
    if contact_dictionary.has_key(key):    
        ### The contact state has already been found
        return contact_dictionary[key]

    else:    
        ### create a new dictionary entry, and record it in the .cindex file
        indexnum = len(contact_dictionary)      
        contact_dictionary[key] = indexnum
        return contact_dictionary[key]


def is_gzipped( filename ):

     fields = splitfields( filename, '.' )
     if len(fields) < 1:
        return False
     else:
         if (fields[ len(fields)-1 ] == 'gz'):
             return True
         else:
             return False



def pairs( free_energies, TargetTemp, statefile, indexfile, OutputPath, PairList ):



    ### Parameters

    kT = K_BOLTZ*TargetTemp
    beta = 1.0/kT
    
    FREE_ENERGY_THRESHOLD = 300.0       # ignore everyting larger than (kcal/mol)

    F_lines_dict = {}
    F_lines = []

    fstates = gzip.open(statefile,'r')
    lines =  fstates.readlines()  
    fstates.close()     
    lines.pop(0)        # get rid of header
    for line in lines:
       fields = splitfields(line)
       F_lines_dict[ int( fields[0] ) ] = float(fields[1])      
       F_lines.append( fields[1] )

    cstates_dict = {}
    findex = gzip.open(indexfile,'r')
    lines =  findex.readlines()
    lines.pop(0)        # get rid of header
    findex.close()
    for line in lines:
        fields = splitfields(line)
        k = int(fields.pop(0))        # get rid of  index
        cstates_dict[ k ] = []
        fields.pop(0) # get rid of  '['
        while len(fields) > 1:
          i = int(fields.pop(0))
          j = int(fields.pop(0))
          cstates_dict[ k ].append( (i,j) )
        
    if DEBUG:
        print 'Number of free energies:',len(F_lines_dict)
        print 'Number of contact states:',len(cstates_dict)



    #####################
    # Calculate the coops

    total = len(F_lines_dict)
    print 'Before filtering,',len(cstates_dict),'contact states'

    ### Get rid of all the infinities and those above threshold from the list
    print 'Removing contact states with F=inf and above threshold',FREE_ENERGY_THRESHOLD,'...'
    i = 0
    for i in F_lines_dict.keys():

        if (F_lines[i] == 'Inf') | (F_lines[i] == 'inf'):
            del F_lines_dict[i]
            del cstates_dict[i]
        if float(F_lines_dict[i]) > FREE_ENERGY_THRESHOLD:
            del F_lines_dict[i]
            del cstates_dict[i]


    print 'After filtering,',len(cstates_dict),'contact states'

    print 'Calculating cooperativities...'


    coop2data = calc_doublecoops(F_lines_dict, cstates_dict, OutputPath, beta)
    
    return coop2data



#################
#
#  class structure for a Dijkstra algorithm vertex (used in ECO calculation)

class Vertex:

    # This is a vertex class for implementing Dijkstra's algorithm
    # a la Cormen et al,   
    
    def __init__(self, res):
    
        self.d = 999    # the shortest path estimate
        self.pi = None  # predecessor node, either a vertex, or None

        self.res = res  # the residue number    

    def __repr__(self):
    
        return str(self.res) 


#############
# Functions #
#############

def cull_contacts(cstates_dict,allcontacts):
 
    
    print 'Culling contact states...'
    
    for key in cstates_dict.keys():
        
        if key % 100 == 0:
            print '\tExamined',key,'of', len(cstates_dict.keys()),'contact states...'
        for c in cstates_dict[key]:
            
            if allcontacts.count(c) == 0:
                allcontacts.append(c)

    if VERBOSE: print '...Done.'
    
    

def calc_doublecoops(F_lines_dict, cstates_dict, OutputPath, beta):
 

    
    ### Cull all contacts
    
    allcontacts = []
    cull_contacts(cstates_dict,allcontacts)
    
    if (DEBUG):
        print 'There are',len(allcontacts),'possible contacts:'
        print allcontacts
        print F_lines_dict
        print cstates_dict
    
    
    outfile = OutputPath + '/coop.txt.gz'
    probfile = OutputPath + '/cprobs.txt.gz'

    print 'Writing two-contact coops to',outfile,'...'    
    fdata = gzip.open(outfile,'w')
    fdata.write( 'c1_i    c1_j    c2_i    c2_j    coop      stab      ECOpair\n')

    
    ### First calculate 
    ### 1) partition function over everything
    ### 2) all single-contact free energies ==> probabilies
    ###    and store them in a dictionary
    
    contactpair_dict = {}       # indexed of 0 through (ncontacts \choose 2)-1
    contact_dict = {}           # indexed of 0 through (ncontacts-1)
    Z = 0.0
    
    contact_lookup = {}         # nested dict to give a quick lookup of contact index
    pair_lookup = {}            # nested dict to give a quick lookup of pair index
    
    ### Here's the trick to make it faster:
    ###
    ### Instead of looping over all contact pairs and going through the list of 
    ### F_states each time,  go through F_states once, and calculate 
    ### everything we need along the way.
    
    
    good_contacts = allcontacts
    good_contacts.sort()
    
    ncontacts = len(good_contacts) 
    ncontact_pairs = ncontacts*(ncontacts-1)/2
    
        
    if DEBUG:
        print 'len(good_contacts)', ncontacts
        print 'good_contacts',good_contacts
        print 'CO_LIMIT = ', CO_LIMIT
     
        print 'Initializing Zcs[] and Zpairs[] lists:'
    
    Zcs = []
    for i in range(0, ncontacts):
        Zcs.append(0.0)
    Zpairs = []         # p(1,1)
    for i in range(0, ncontact_pairs):
        Zpairs.append(0.0)

    if DEBUG:
        print 'Initializing contact_lookup dict...'
    for i in range(0, ncontacts):
      contact_lookup[ good_contacts[i] ] = i 
      contact_dict[ i ] = good_contacts[i] 
        
    
    if DEBUG:
        print 'Initializing pair_lookup dicts...'
    pairindex = 0
    for i in range(0, ncontacts-1):
      for j in range(i+1, ncontacts):
        
        if pair_lookup.has_key( good_contacts[i] ):     
            pair_lookup[good_contacts[i]][good_contacts[j]] = pairindex
        else:   
            pair_lookup[good_contacts[i]] = {}
            pair_lookup[good_contacts[i]][good_contacts[j]] = pairindex

        pairindex = pairindex + 1
        
    # print 'pair_lookup',pair_lookup   
        
        
        
    if DEBUG:
        print 'Compiling pair contacts Zcs[] and Zpairs[] for each contact state...'
    keycount = 0
    nkeys = len(F_lines_dict.keys())
    for key in F_lines_dict.keys(): 

        keycount = keycount + 1
        if (keycount % 500) == 0:
           print keycount,'of', nkeys,'states...'
        
         
        thisfreeenergy = F_lines_dict[key]
        this_cstate = cstates_dict[key]
        
        if (thisfreeenergy != 'Inf') & (thisfreeenergy != 'inf') & (thisfreeenergy != None) :   

            thisZ = exp(-1.0*beta*thisfreeenergy)
            Z = Z + thisZ
            
            if DEBUG:
                print 'Examining state',key, this_cstate, thisZ
                            
            for i in range(0, len(this_cstate)):
            
              c1 = this_cstate[i]
              if good_contacts.count(c1) > 0:         
                  Zcs[ contact_lookup[c1] ] = Zcs[ contact_lookup[c1] ] + thisZ
                  if DEBUG:
                      print '\tAdding',thisZ,'to Zc[',c1,']'
              
                  for j in range(i+1, len(this_cstate)):
              
                    c2 = this_cstate[j]
                
                    if good_contacts.count(c2) > 0:
                        pairindex = pair_lookup[c1][c2]
                        # print c1,c2,pairindex

                        Zpairs[pairindex] = Zpairs[pairindex] + thisZ
          
            
    
    # print Zcs
    # print Zpairs
    
    
    print 'Writing contact probabilities to file',probfile,'...'
    fprob = gzip.open(probfile,'w')     
    fprob.write( 'ci      cj      prob      co\n' )             
    for i in range(0, ncontacts):    
        thiscontact = contact_dict[i]
        thisprob = Zcs[i]/Z
        co = thiscontact[1] - thiscontact[0]
        if co > CO_LIMIT:
            tmpstr = str(thiscontact[0])
            while len(tmpstr) < 8:
                tmpstr=tmpstr+' '
            tmpstr = tmpstr+str(thiscontact[1])
            while len(tmpstr) < 16:
                tmpstr=tmpstr+' '               
            fprob.write( tmpstr )
            tmpstr = '%1.10f'%thisprob
            fprob.write( tmpstr[0:7]+'   ' )
            fprob.write(str(co) + '\n' )
                
    fprob.close()                

    coop2data = zeros( (len(Zpairs), 7), float)
    
    for c1 in pair_lookup.keys():

      c1_index = contact_lookup[c1]
      
      for c2 in pair_lookup[c1].keys():

          pairindex = pair_lookup[c1][c2]
          c2_index = contact_lookup[c2]
          
          ### Calculate the relative entropy of joint vs. prod of indep. probs

          D = 0.0        # the relative entropy

          p_c1 = Zcs[c1_index]/Z
          p_c2 = Zcs[c2_index]/Z
          p_notc1 = 1 - p_c1
          p_notc2 = 1 - p_c2

          p11 = Zpairs[pairindex]/Z
          p10 = p_c1 - p11
          p01 = p_c2 - p11
          p00 = 1.0 - p11 - p10 - p01

          ### notc1, notc2
          if p00 > 0.0:
            if (p_notc1*p_notc2) > 0.0:
              Dterm = (p00) * math.log( (p00/( p_notc1 * p_notc2)) ) * INVLOG2
              if DEBUG:
                  print 'D[notc1][notc2] =',Dterm
              D = D + Dterm

          ### c1, notc2
          if p10 > 0.0:
            if (p_c1*p_notc2) > 0.0:
              Dterm =  (p10 ) * math.log( p10/( p_c1 * p_notc2) ) * INVLOG2
              if DEBUG:
                  print 'D[c1][notc2] =',Dterm
              D = D + Dterm

          ### notc1, c2
          if p01 > 0.0:
            if (p_notc1*p_c2) > 0.0:
              Dterm = (p01) * math.log( p01/( p_notc1 * p_c2) ) * INVLOG2
              if DEBUG:
                  print 'D[notc1][c2] =',Dterm
              D = D + Dterm

          ### c1, c2
          if p11 > 0.0:
            if (p_c1*p_c2) > 0.0:
              Dterm = (p11) * math.log( p11/( p_c1 * p_c2) ) * INVLOG2
              if DEBUG:
                  print 'D[c1][c2] =',Dterm
              D = D + Dterm


          stability = p11
          if stability < p_c1*p_c2:     # anti-cooperative?
            D = -1.0*D

          ecos = []
          ecos.append(eco_newcontact([c1],c2))
          ecos.append(eco_newcontact([c2],c1))
          if DEBUG:
                  print 'ECO',c1,'+',c2,ecos[0]
                  print 'ECO',c2,'+',c1,ecos[1]
          mineco = min(ecos)
          if DEBUG:
                  print 'minECO',mineco
          
          
          ### write output in 8-char wide fields (NO TABS, for scott's sake)
          t = str(c1[0])
          while len(t) < 8:
              t=t+' '
          t = t + str(c1[1])
          while len(t) < 16:
              t=t+' '
          t = t + str(c2[0])
          while len(t) < 24:
              t=t+' '
          t = t + str(c2[1])
          while len(t) < 32:
              t=t+' '
          fdata.write( t )
          
          t = '%2.10f'%D
          fdata.write( t[0:8]+'  ' )
            
          t = '%2.10f'%stability
          fdata.write( t[0:8]+'  ' )
          
          t = str(mineco)
          while len(t) < 32:
              t=t+' '
          fdata.write( t + '\n')
          
                  

          if DEBUG:
              print c1,c2,'\tD:',D
          
          coop2data[pairindex,0] = c1[0]   
          coop2data[pairindex,1] = c1[1]   
          coop2data[pairindex,2] = c2[0]   
          coop2data[pairindex,3] = c2[1]   
          coop2data[pairindex,4] = D   
          coop2data[pairindex,5] = stability   
          coop2data[pairindex,6] = mineco

               

    if VERBOSE: print '...Done.'
    fdata.close()
    
    return coop2data



def eco_newcontact(knowncontacts,newcontact):


        # calculates the ECO=tau of a new contact (k,l) with respect to existing contacts 
        # using Dijkstra's algorithm
        #
        # Procedure: 
        #
        #     1) create the eco_matrix with eco=1 for known contacts, and CO for all others  
        #     2) implement Dijkstra searching    

        if len(knowncontacts)==0:    # if there's just one contact return the CO
            k = newcontact[0]
            l = newcontact[1]
            return abs(k-l)

        if len(knowncontacts)>0:

          # First -- initialize eco_matrix(i,j)
          # find the maximum residue number for matrix size
          
          maxres = 0
          for contact in knowncontacts:
              if contact[1] > maxres:    # last tuple is always largest 
                  maxres = contact[1]
          if newcontact[1] > maxres:    # last tuple is always largest 
              maxres = newcontact[1]

          eco_matrix = []    # initialize eco_matrix with COs
          for i in range(0,maxres+1):
              eco_matrix.append([])
              for j in range(0,maxres+1):    # symmetirc to avoid confusion
                eco_matrix[i].append(abs(i-j))    # for initialization

          for contact in knowncontacts:
              eco_matrix[contact[0]][contact[1]] = 1    # these are in contact
              
          
          ### Construct the Dijkstra graph of all residues
          allres = []
          allcontacts = []
          for contact in knowncontacts:
              allcontacts.append(contact)
          allcontacts.append(newcontact)
          
          for contact in allcontacts:         
              if allres.count(contact[0])==0:
                  allres.append(contact[0])             
              if allres.count(contact[1])==0:
                  allres.append(contact[1])
              if allres.count(contact[0])==0:
                  allres.append(contact[0])             
              if allres.count(contact[1])==0:
                  allres.append(contact[1])
          allres.sort()         # sort in place, for easier debugging later     
          # print 'allres',allres  
            
          graph = []
          for res in allres:
              graph.append(Vertex(res))
          # print 'graph',graph    
                        
         
          ### Find the vertex objects that correspond to the newcontact (k,l)
          for i in range(0,len(graph)):
              if graph[i].res == newcontact[0]:
                  k = graph[i]
              if graph[i].res == newcontact[1]:
                  l = graph[i]
          # print 'k,l',k,l
          ### Find the vertex objects that correspond to the newcontact (k,l)
          eco = dijkstra_eco(graph,eco_matrix,k,l)
          
          return eco


def dijkstra_eco(graph,eco_matrix,source,end):

    ### Returns ECO of (source,end) using
    ### Dijkstra's algorithm to find the shortest paths to all vertices
    ### given a starting source node

    initialize_single_source(graph,source)

    s = []  # the set of vertices for which we have already computed
                # the final shortest path.

    queue = []

    # add all vertices to the queue 
    for v in graph:
        queue.append(v)


    while len(queue) > 0:           ## the first time through,
                                        ## the source will be u

        u = extract_min(queue)      ### pop off the         
        s.append(u)

        for v in graph:

            ### For the ECO problem, a vertex is adjacent to all other vertices
            if v != u:

                relax(u,v,eco_matrix)

    return end.d


def initialize_single_source(graph,source):

    ### a graph is just a list of vertices

    for v in graph:

        v.d = 999
        v.pi = None

    source.d = 0    

def extract_min(queue):

    ### pops off a vertex v from the queue with the smallest v.d

    bestd = 999
    for i in range(0,len(queue)):

        if queue[i].d < bestd:
            bestd = queue[i].d
            besti = i

    ### pop off the best d and return it
    return queue.pop(besti)


def relax(u,v,eco_matrix):

    ### re-adjusts the upper bound on the known shortest distance d

    if v.d > ( u.d + eco_matrix[u.res][v.res] ):   
        v.d = u.d + eco_matrix[u.res][v.res]
        v.pi = u





if __name__ == "__main__":

  if len(sys.argv) < 5:
      print usage
    
  else:

      statefile = sys.argv[1]  
      cindexfile = sys.argv[2]  
      ckeyfile = sys.argv[3]  
      outfile = sys.argv[4]
  
      SELECT_FLAG = 0
      select_residues = []
      select_strings = []
      absorb_contacts = []
      contact_strings = []


      CO_LIMIT = 0
      for i in range(0,len(sys.argv)):
          if sys.argv[i] == '-CO':
              CO_LIMIT = int(sys.argv[(i+1)])
              sys.argv.pop(i)
              sys.argv.pop(i)
              break


      print 'sys.argv', sys.argv
      for i in range(0,len(sys.argv)):
          if sys.argv[i] == '-select':
              SELECT_FLAG = 1
              for j in range(0,len(sys.argv[(i+1):])):
                  select_strings.append( sys.argv.pop() )
              sys.argv.pop()   # get rid of -select
              break

          print 'select_strings',select_strings

          if sys.argv[i] == '-absorb':
              contact_strings = sys.argv[(i+1):(len(sys.argv)+1)]
              break

      print 'contact_strings',contact_strings
      if len(contact_strings) % 2 != 0:
         print 'Error:  An even number of residues needed for contact pairs. Exiting.'
         sys.exit(0)

      i=0
      while i < len(contact_strings):
          absorb_contacts.append( ( int(contact_strings[i]),  int(contact_strings[i+1]) ) )
          i = i + 2

      i=0
      while i < len(select_strings):
          select_residues.append( int(select_strings[i]) )
          i = i + 1

      select_residues.sort()
      print 'select_residues:', select_residues    

      config = Config('for pairs')
      pairs( config, statefile, cindexfile, ckeyfile, outfile, CO_LIMIT, absorb_contacts, SELECT_FLAG, select_residues)
      

      
def punch( coop2data, OutputPath, PairList, CritCutProb,  CritCutCoop ):
        
    print 'Calculating punch scores for all contacts...'
    
    outfile = OutputPath + '/punch.txt'
    fout = file(outfile,'w')
        
    nres = int(max(  max(coop2data[:,1]), max(coop2data[:,3]) )) + 1
    # print 'nres:',nres

    c_degree = zeros( (nres,nres), float ) 
    s_degree = zeros( (nres,nres), float ) 


    coop2dims = shape(coop2data)
    # print 'coop2dims',coop2dims
    npairs = coop2dims[0]

    # Recall, an empty pairlist means examine all conatcts
    if len(PairList) == 0:
      for i in range(0,npairs):
        c1 = (int(coop2data[i,0]),int(coop2data[i,1]))
        c2 = (int(coop2data[i,2]),int(coop2data[i,3]))
        if PairList.count( c1 ) == 0:
            PairList.append( c1 )
        if PairList.count( c2 ) == 0:
            PairList.append( c2 )

    
    punchlist = []
    PunchStabVals = []
    PunchCoopVals = []
    
    # print PairList
    for c in PairList:

        I_connected = ( (coop2data[:,0] == c[0])&(coop2data[:,1] == c[1]) ) | ( (coop2data[:,2] == c[0])&(coop2data[:,3] == c[1]) )
        # print 'I_connected',I_connected
        
        coop2data_connected = coop2data[I_connected]
        # print 'coop2data_connected',coop2data_connected

        c_sorted = sort(coop2data_connected[:,4])
        c_sorted = c_sorted[::-1]
        I_degree = c_sorted > CritCutCoop
        c_punch = float(len(c_sorted[I_degree]))/float(len(PairList))
        c_degree[c[0],c[1]] = c_punch
        
        s_sorted = sort(coop2data_connected[:,5])
        s_sorted = s_sorted[::-1]
        I_degree = s_sorted > CritCutProb
        s_punch = float(len(s_sorted[I_degree]))/float(len(PairList))
        s_degree[c[0],c[1]] = s_punch
        
        punchlist.append( [s_punch, c_punch, c] )
        PunchStabVals.append(s_punch)
        PunchCoopVals.append(c_punch)
        # print c, punchlist


        
    PunchStabTot = sum(sum(s_degree))    
    PunchCoopTot = sum(sum(c_degree))
    
    fout.write('PUNCH STABILITY TOTAL:\n')
    fout.write(str(PunchStabTot)+'\n\n')
    fout.write('PUNCH COOPERATIVITY TOTAL:\n')
    fout.write(str(PunchCoopTot)+'\n\n')
    
    punchlist.sort()
    punchlist.reverse()

    if DEBUG:
        print 'punchlist',punchlist
        
    fout.write('CONTACT PUNCH VALUES:\n')
    fout.write('contactname stab_punch  coop_punch\n')
    for entry in punchlist:
    
    
        c = entry[2]
        cstr = ("%d,%d" % c).ljust(12)    
        fout.write(cstr)

        stab_punch = entry[0]
        coop_punch = entry[1]
                
        tmpstr = '%2.10f'%stab_punch
        fout.write(tmpstr[0:9]+'   ')
    
        tmpstr = '%2.10f'%coop_punch
        fout.write(tmpstr[0:9]+'\n')   
    
    
    if VERBOSE: print '...Done'
    return [PunchStabTot, PunchCoopTot, PunchStabVals, PunchCoopVals]

