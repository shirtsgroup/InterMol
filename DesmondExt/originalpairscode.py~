          elif match.group('pairs'): #GromacsTopology didn't fix this, so fix this when changes are made
            if verbose:
              print "Parsing [ pairs]..."
            print "pairs at line %d"%i
	    i+=6 #ASSUMING ALL COLUMNS OF PAIRS ARE THE SAME
	    funct1_lj = []
	    funct1_cl = []
	    funct1_both = []
	    funct2 = [] 
	    ljch=None
	    clomb=None
	    lj = False
	    cl = False
            tempcnt = 0
	    #j = i
	    while not re.match(r'\s*[:::]',lines[i]): #organzies pair data into 4 types above
	      split = lines[i].split()
	      if int(split[0]) == 1:
	        if re.match(split[3], "LJ"):
		  ljch = float(split[4])
		  lj = True
		  funct1_lj.append([int(split[1]),int(split[2])])
		else:
		  clomb = float(split[4])
		  cl = True
		  funct1_cl.append([int(split[1]),int(split[2])])
	      elif lj:
	        if not float(split[4]) == ljch and re.match(split[3], "LJ"):
		  funct2.append([int(split[1]),int(split[2]), split[3],float(split[4])])
		elif re.match(split[3], "Coulomb"):
	          if not clomb:
                    clomb = float(split[4])
		  cl = True
		  lj = False 
	      elif cl:
	        if not float(split[4]) == clomb and re.match(split[3], "Coulomb"):
		  funct2.append([int(split[1]),int(split[2]), split[3],float(split[4])])
		elif re.match(split[3], "LJ"):
	          if not ljch:
                    ljch = float(split[4])
		  cl = False
		  lj = True
	      if tempcnt == len(funct2) and not int(split[0]) == 1:
	        if lj and [int(split[1]), int(split[2])] in funct1_cl:
	  	  funct1_cl.remove([int(split[1]), int(split[2])])
		  funct1_both.append([int(split[1]), int(split[2])])
	        elif cl and [int(split[1]), int(split[2])] in funct1_lj:
	          funct1_lj.remove([int(split[1]), int(split[2])])
		  funct1_both.append([int(split[1]), int(split[2])])
	        else:
		  if lj:
		    funct1_lj.append([int(split[1]),int(split[2])])
		  else:
		    funct1_cl.append([int(split[1]),int(split[2])])
	      else:
		if int(split[0]) > 1:
	          tempcnt+=1		
              i+=1
	    #j = 0
	    j = 0
	    for a in funct2: #PUT ALL PAIRS WITH FUNCT 2 INT LJ2
	      newPairForce = AbstractPair(a[0], a[1], 2) 
              currentMoleculeType.pairForceSet.add(newPairForce)
	      System._sys._forces.add(newPairForce)
	      
	      if System._sys._combinationRule == 1:
                newPairType = LJ2PairCR1Type(atomlist[a[0]-1].atomName, #atom 1 and index
                                atomlist[a[1]-1].atomName,              #atom 2 and index
                                2,                                     #type
				1,                                     #fudgeqq
				float(ljch) * units.elementary_charge, #UNITS CORRECT???
                                float(clomb) * units.elementary_charge, #UNITS CORRECT???
                                float(0) * units.kilojoules_per_mole * units.nanometers**(6),  #COME BACK
                                float(0) * units.kilojoules_per_mole * units.nanometers**(12)) #COME BACK
	      elif System._sys._combinationRule == (2 or 3):
                newPairType = LJ2PairCR23Type(atomlist[a[0]-1].atomName,
                                atomlist[a[1]-1].atomName,
                                2,
				1,
				float(ljch) * units.elementary_charge, #UNITS CORRECT???
                                float(clomb) * units.elementary_charge, #UNITS CORRECT???
                                float(0) * units.nanometers,          #COME BACK
                                float(0) * units.kilojoules_per_mole) #COME BACK
	      self.pairtypes.add(newPairType)
	      
	    
	    for a in funct1_lj:     #PUT ALL PAIRS WITH ONLY LJ INTO LJ1
	      newPairForce = AbstractPair(a[0], a[1], 1) 
              currentMoleculeType.pairForceSet.add(newPairForce)
	      System._sys._forces.add(newPairForce)
	      
	      if System._sys._combinationRule == 1:
                newPairType = LJ1PairCR1Type(atomlist[a[0]-1].atomName,  #atom 1 and index
                                atomlist[a[1]-1].atomName,              #atom 2 and index
                                1,                                     #type
                                float(0) * units.kilojoules_per_mole * units.nanometers**(6),  #COME BACK
                                float(0) * units.kilojoules_per_mole * units.nanometers**(12)) #COME BACK
	      elif System._sys._combinationRule == (2 or 3):
                newPairType = LJ1PairCR23Type(atomlist[a[0]-1].atomName,
                                atomlist[a[1]-1].atomName,
				1,
                                float(0) * units.nanometers,          #COME BACK
                                float(0) * units.kilojoules_per_mole) #COME BACK
	      self.pairtypes.add(newPairType)
	   
	    for a in funct1_cl: #PUT ALL PAIRS WITH ONLY CL INT LJ1
	      newPairForce = AbstractPair(a[0], a[1], 1) 
              currentMoleculeType.pairForceSet.add(newPairForce)
	      System._sys._forces.add(newPairForce)
	      
	      if System._sys._combinationRule == 1:
                newPairType = LJ1PairCR1Type(atomlist[a[0]-1].atomName,  #atom 1 and index
                                atomlist[a[1]-1].atomName,               #atom 2 and index
                                1,                                     #type
                                float(0) * units.kilojoules_per_mole * units.nanometers**(6),  #COME BACK
                                float(0) * units.kilojoules_per_mole * units.nanometers**(12)) #COME BACK
	      elif System._sys._combinationRule == (2 or 3):
                newPairType = LJ1PairCR23Type(atomlist[a[0]-1].atomName,
                                atomlist[a[1]-1].atomName,
				1,
                                float(0) * units.nanometers,          #COME BACK
                                float(0) * units.kilojoules_per_mole) #COME BACK
	      self.pairtypes.add(newPairType)
	     
	    for a in funct1_both: #PUT ALL PAIRS WITH BOTH LJ AND COULOMB INTO NB 
	      newPairForce = AbstractPair(a[0], a[1], 1) 
              currentMoleculeType.pairForceSet.add(newPairForce)
	      System._sys._forces.add(newPairForce)
	      
	      if System._sys._combinationRule == 1:
                newPairType = LJNBPairCR1Type(atomlist[a[0]-1].atomName, #atom 1 and index
                                atomlist[a[1]-1].atomName,               #atom 2 and index
                                2,                                     #type
				float(ljch) * units.elementary_charge, #UNITS CORRECT???
                                float(clomb) * units.elementary_charge, #UNITS CORRECT???
                                float(0) * units.kilojoules_per_mole * units.nanometers**(6),  #COME BACK
                                float(0) * units.kilojoules_per_mole * units.nanometers**(12)) #COME BACK
	      elif System._sys._combinationRule == (2 or 3):
                newPairType = LJNBPairCR23Type(atomlist[a[0]-1].atomName,
                                atomlist[a[1]-1].atomName,
                                2,
				float(ljch) * units.elementary_charge, #UNITS CORRECT???
                                float(clomb) * units.elementary_charge, #UNITS CORRECT???
                                float(0) * units.nanometers,          #COME BACK
                                float(0) * units.kilojoules_per_mole) #COME BACK
	      self.pairtypes.add(newPairType)
	      #if j < 20000:
              #  print "Atom %s and Atom %s in Row %d"%(atomlist[a[0]-1].atomName,atomlist[a[1]-1].atomName,j)
	      #j+=1 