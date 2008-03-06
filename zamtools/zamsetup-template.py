#!/usr/bin/env python

#LAST MODIFIED: 05-15-07


#NOTE: All of the variables below represent the default values
#for possible changes.  DELETE any variables you are not changing
#as these may override settings or reset data that have evolved
#in the zam simulation.


#========================= BASIC CHANGES =========================

zd.PreBuildString = ""                   #string run in tleap before system ('sys') is built [string]
zd.PostBuildString = ""                  #string run in tleap after system ('sys') is built [string]
zd.UpdateFrags = False                   #indicates whether or not to automatically create new fragments [boolean]
zd.UpdateRest = False                    #whether or not to update master restraints automatically, not currently used [boolean]

ConserveDisk = False   #True will delete large files to save space
ShowRest = False       #whether or not to show proposed restraints
ExcessClients = 1      #excess number of socket clients to keep



#====== MODULE CHANGES ======


#==== zamdata.py variables ====
zamdata.PairMinCO = 3        #Minimum CO for pair analyses
zamdata.PairMaxCO = 1000     #Maximum CO for pair analyses
zamdata.PairMinECO = 3       #Minimum ECO for considering a pair for analysis
zamdata.AllPairsDflt = True  #analyze all (true) or just hydrophobic (false) pairs  


#==== zamchoice.py variables ====
zamchoice.NResGap = 1             #stride in num of residues between initial 8mer frags
zamchoice.MaxFragsDflt = 10       #maximum number of new fragments to create in one round after 16mers
zamchoice.MaxChoicesDflt = 100    #maximum number of assembly choices to report
zamchoice.NPredictions = 10       #maximum number of full chain predictions to make
zamchoice.NResLoop = 4            #won't penalize loops of this length or less
zamchoice.AssemResPenalty = 5.0   #penalty per res per num of previous assemblies in score 
zamchoice.COWeight = 10.75        #weight multiplied by contact order and added to score


#==== zamrex.py variables ====
#global replica exchange run variables
zamrex.NRepIncr = 5       #will only have integer multiples of this many replicas;
                   #a negative value uses the current number of processors
                   #with a minimum of abs(NRepIncr)
zamrex.MinTemp = 270.0    #minimum replica temperature 
zamrex.MaxTemp = 600.0    #maximum replica temperature
zamrex.SwapRest = True    #swap restraints when they are sorted?

#cycle duration; first number is fragment size, second is cycle time in ps
#for all fragments of that size or less. fragments must be in order
#smallest to largest.  
zamrex.CycleTimes = [(8, 20.), (12, 10.), (16, 5.), (24, 2.), (30, 1.)]

#rigidifying restraint params
zamrex.RestRigidFConst = 0.1

#secondary structure restraint params
#these are (Phi, Psi, PhiTol, PsiTol)
zamrex.RestSSFConst = 0.5
zamrex.RestSSParams = {"H":(-58., -47., 10., 10.),
                "E":(-130., 125., 40., 40.)}

#contact restraint params; the distances here are absolute
zamrex.RestContactParams = {"DIST1":0., "DIST2":0., "DIST3":7.5, "DIST4":8.0,
    "FCONST2":0.0, "FCONST3":0.5}
zamrex.RestContAtom = "*"  #use residue centroid

#positional restraints
zamrex.PosRestFConst = 1.0

#ANCHORING RESTRAINTS FOR MULTIPLE CHAINS:
#this is the distance added to Rsurf1+Rsurf2 at which a restraint starts to
#hold chains near each other so they don't fly off into space
zamrex.ChainHoldDist = 10.
zamrex.ChainHoldFConst = 1.0

#UMBRELLA SAMPLING RESTRAINTS
#umbrella sampling force constant
zamrex.UmbrFConst = 0.5
#width of umbrella bins
zamrex.UmbrDelta = 5.0
#umbrella atom
zamrex.UmbrAtom = "*"  #centroid

#temperature optimization
zamrex.TempOptTimes = []    #times in ps at which temperatures are optimized

#minimum fractional restraint strength when scaling
zamrex.ScaleRestMin = 0.


#==== zamrun.py variables ====
zamrun.MaxInitConf = 50      #maximum number of initial conformations
zamrun.MaxInitAssem = 10     #maximum number of pairwise assemblies to attempt
zamrun.RexTimeBase = 5000.   #baseline time in ps to REX
zamrun.RexTimePerAdd = 0.    #addl time in ps to REX per added residue
zamrun.RexTimePerDOF = 0.    #addl time in ps to REX per degree of freedom
zamrun.RexTimeMin = 4000.    #minimum total time to REX for any fragment
zamrun.RexTimeSwapRest = 0.  #initial time to anneal contacts for wait fragments
zamrun.AnalTime = 1000.      #time in ps to use for analysis
zamrun.InitConfTimeout = None #timeout in seconds for the initial configuration making, or None


#==== zamanal.py variables ====

#analysis variables
zamanal.DistMethod = 1      #contact distances calc between 0 = CA, 1 = CB, 2 = residue centroid
zamanal.PmfNBin = 11        #num of pmf bins
zamanal.PmfBinMin = 4.      #minimum pmf distance
zamanal.PmfBinMax = 15.     #maximum pmf distance
zamanal.ContactCut = 8.0    #distance defining a contact cutoff
zamanal.TargetTemp = 270.   #target temperature for calculations
zamanal.DeltaTemp = 100.    #maximum temperature range to consider in WHAM calculations
zamanal.DeltaTempUmbr = 200.#maximum temperature range to consider in WHAM when umbrella sampling
zamanal.RemoveRest = True   #remove restraints using wham?
zamanal.ChainPmfNBin = 46   #number of bins in chain pmf
zamanal.ChainPmfMin = 4.    #minimum chain pmf distance
zamanal.ChainPmfMax = 50    #maximum chain pmf distance
zamanal.UmbrNBin = 50       #num of pmf bins for umbrella pmfs
zamanal.UmbrWindow = 10.    #num of angstroms on each side of umbrella range to do min, max for pmf
zamanal.UmbrDistMethod = 1  #umbrella distances calc between 0 = CA, 1 = CB, 2 = residue centroid

#clustering parameters
zamanal.MaxCluster = 10     #maximum number of cluster configurations
zamanal.ClustRmsdTol = 2.0  #tolerance rmsd between clusters
zamanal.ClustMaxIter = 10   #maximum number of iterations

#punch analysis
zamanal.PunchCritProb = 0.5 #probability threshold for punch calculation
zamanal.PunchCritCoop = 0.4 #cooperativity threshold for punch calculation

#score parameters
zamanal.PmfCutDist = 9.0     #distance in angstroms for taking pmf difference
zamanal.ScoreSlope = -14.0   #slope for critical scores in (Cut,PMFDiff) plane
zamanal.ScoreInt = 11.0      #intercept for critical scores in (Cut,PMFDiff) plane

#bestn analysis
zamanal.BestnMinECO = 3      #minimum ECO for generating bestn pairs
zamanal.BestnN = 3           #number of contact pairs to choose

#energy analysis
zamanal.EneTerms = ["TEMP","ETOT","EPOT","EREST","EVDW","EEL","ESURF"]
                     #terms to average from the energy function

#score analysis
zamanal.ScoreNStride = 20    #stride in frames for calculating various scores

#timeouts
zamanal.AnalTimeout = None #timeout in seconds for the analysis routines
zamanal.SnapTimeout = None #timeout in seconds for the snapshot routines


#==== zamcontact.py variables ====
zamcontact.UseCO = False
zamcontact.MinPropRest = 1
zamcontact.MaxPropRest = 3
#these are coefficients for (prior, cprob, pmfscore, punchstab, punchcoop, meso)
zamcontact.LogitParams = [(8,  0, -2.5916,  2.9295,  0.0000, -0.0820,  0.0000, -0.0787),
                          (12, 0, -2.3686,  2.5255,  0.0000, -0.0373,  0.0000,  0.1154),
                          (16, 0, -1.9864,  2.7035,  0.0107, -0.0262,  0.0000,  0.0306)]
#these add priors and correct for correlations among different measurements;
#these are coefficients for (maxCO, b0, b1, b2, b3, b4)
zamcontact.CorrelParams = [(4, 0.9171, 0.0337, 0.9892, 0.0000, 1.0000, 0.0000),
                           (0, 0.9171, 0.0337, 0.9892, 0.0000, 1.0000, 0.0000)]

