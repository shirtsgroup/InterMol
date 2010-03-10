from mmtools.gromacstools.devel.Scripts.extractNumberSolventMoleculesFromTopology import *
from mmtools.gromacstools.devel.Scripts.IonEnvironmentCalculation import *

def calculateCounterionsToAdd(systemSetup, topfile, totalSystemCharge):
    """Calculate the number of counterions needed to both
    1) neutralize charge and 2) approximate the desired concentration
    
    REQUIRED
    
    systemSetup     - a SystemSetup() object containing the salt conditions for the system
    topfile         - the filename of the topology file for the system (assumed to already be solvated)
    
    RETURNS
    
    np                   - the number of water molecules than should be replaced by positive ions 
    nn                   - the number of water molecules than should be replaced by negative ions 
    numberWaterMolecules - the number of water molecules currently in the simulation box
    
    
    NOTES
    
    
    Here's a description of the method we use to do this:
    
    # Our philosophy here will adhere to simple rules:
    #     1) always neutralize the charge in the box
    #     2) always decrease the total numbers of ions, if possible
    #     3) there will be divalent anions, only -1 charges.
    
    ######################################################################################## 
    # Let f be the fraction of salt molecules in a solution of total solvent molecules N.  #
    # When the ions dissociate, there are a total of M*N*f (cations+anions), where         #
    #                                                                                      #
    #     M = (systemSetup.positiveStoichiometry + systemSetup.negativeStoichiometry)      #
    #                                                                                      #
    # When we add cations an anions with the genion* program, each *ion* displaces a       #
    # solvent molecule! How many total salt molecules x should we introduce to maintain    #
    # the correct concentration f?  x must satisfy the equation:                           #
    #                                                                                      #
    #     x/(N-Mx) = f                                                                     #
    #                                                                                      #
    # which has the solution:                                                              #
    #                                                                                      #
    #     x = fN/(1+Mf)                                                                    #
    #                                                                                      #
    ########################################################################################
    """

    # Estimate to the nearest integer charge
    totalSystemCharge = int(totalSystemCharge)
    
    # Extract the number of water molecules from the topology file
    numberWaterMolecules = extractNumberSolventMoleculesFromTopology( topfile )
    
    # create an IonEnvironmentCalculation object for doing ion calculations  
    ioncalc = IonEnvironmentCalculation(systemSetup)
    
    # Calculate the number of ions that should be in solution based on the concentration
    numberFraction = ioncalc.molarityToNumberFraction( systemSetup.saltconc )
    print 'Ratio of salt molecules to solvent molecules:', numberFraction
    
    M = (systemSetup.positiveStoichiometry + systemSetup.negativeStoichiometry)
    numSaltMolecules = int(numberFraction*float(numberWaterMolecules)/(1.+float(M)*numberFraction))
    print '\nNumber of salt molecules to add:', numSaltMolecules
     
    # Start with the number of counterions predicted from salt concentration
    np = systemSetup.positiveStoichiometry*numSaltMolecules
    nn = systemSetup.negativeStoichiometry*numSaltMolecules
    print "\nFinding the best numbers of addded positive and negative ions..."
    
    
    # Incrementally remove charge from the system until we have acheived neutralitry
    totalChargeIsNeutral=False 
    totalIonCharge = (nn*systemSetup.negativeIonCharge + np*systemSetup.positiveIonCharge)                   
    while not totalChargeIsNeutral:
	
	print 'totalIonCharge =',totalIonCharge, '(%d*%d + %d*%d)'%(nn,systemSetup.negativeIonCharge,np,systemSetup.positiveIonCharge),'totalSystemCharge',totalSystemCharge
	
	# If there is too much positive charge, make it more negative
	if (totalIonCharge + totalSystemCharge) > 0:
	    # are there enough pos ions to remove?
	    if np > 0:
		np=np-1
	    else:
		nn=nn+1
	# If there is too much negative charge, make it more positive
	elif (totalIonCharge + totalSystemCharge) < 0:
	    # are there enough neg ions to remove?
	    if nn > 0:
		nn=nn-1
	    else:
		np=np+1
	else:
	    pass
	
	totalIonCharge = (nn*systemSetup.negativeIonCharge + np*systemSetup.positiveIonCharge)    
	
	if (totalIonCharge + totalSystemCharge) == 0:
	    totalChargeIsNeutral=True
    
    print '...Done.'
    print 'Total charge of system before added ions:', totalSystemCharge
    print 'Total charge of added ions:', totalIonCharge
    print '    %d %s ions and %d %s ions needed.'%(nn, systemSetup.negativeIonName, np, systemSetup.positiveIonName)
    
    return [np, nn, numberWaterMolecules]
