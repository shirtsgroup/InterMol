#!/usr/bin/env python

from string import *
from numpy import *

import sys

DEBUG = 0       # turn flag on for detailed statements
                # about what the algorithm is doing....

description = """
There are some things you need to know about bestn.py.   The main workhorse of the bestn.py is an object called ContactGroup.  Here are its methods:

class ContactGroup:
    # An object for calculating and storing the best-n punch hits   
    def __init__( self,  pdict, eco, restrained_contacts, fancy_bounds, restonly):
    def __repr__( self ):
    def too_close( self, ci, cj ):       
    def explore( self, nlimit, selectmax):

explore() is the routine that actually searches for the best N contacts.  

* If restonly == True, then the search will only examine contacts that have ECOs far enough away from the specified restraint contacts.
* If restonly == False, then the search will also try to find sets of contacts that also have ECOs far enough away from *each other*.
-------------
* If fancy_bounds == True, then it will use a "branch-and-bound" method to find the single best set of N contacts
* If fancy_bounds == False, it will try to exhaustively search all possible N-contact selections.

The problem with (fancy_bounds == False) is that the enumeration can, in some cases, take a really long time.  Luckily, I have a fix.  The search starts with the "greedy" strategy that you were using before--picking the highest-scoring contact first, then the next-highest (if it's not too close to the other contacts or restraints), and so on.  
This means that unless the search problem is incredibly diabolical (it's not), most of the time spend enumerating all selections is useless wheel-spinning.  To remedy this I have provided a variable

    self.MAX_ITERS = 100000

that will cut the exhaustive search short, in the case that it starts spinning its wheels.  100000 iterations seems to be a good number (about 20 seconds max).  Your results may vary :)

The other issue I must mention has to do with when (fancy_bounds == True).  In this case, the search is very greedy, so the list of best N-contact groups (best_contact_groups[])  returned by bestn() is not really the best list -- it is only the best groups encountered by the search process.

Here's an example:  Say you do a fragment simulation, with some constraints, and it looks like it wants to be a helix.  For the next round of simulation, you want to put up to three more restraint contacts on based on the punch scores (let's say).  So, you want to generate a list of sets of three or less non-competing (ECO>3) 

You would run bestn() with

    pdict = {<contact>:<punchscore>, ... }
    n = 3
    eco = 3
    restrained_contacts = [ <your restraints> ]
    fancy_bounds == False       (to get a true exhaustive sampling of all possible contact groups)
    restonly == False   (you want to choose a set of contacts that are all far enough away from each other)        

"""

########################
### Classes

class ContactGroup:

    # An object for calculating and storing the best-n punch hits
    
    def __init__( self,  pdict, eco, restrained_contacts, fancy_bounds, restonly):
    
        self.pdict = pdict 
            
        self.contacts = []
        self.weights = []       
        self.select = []
        self.incompatible = {}
        self.restrained_contacts = restrained_contacts
        
        self.incompatible_cutoff = eco
        
        self.selectmax = 0
        self.level = 0
        self.best_score = 0.0
        self.current_score = 0.0
        
        self.best_select = []
        self.best_contacts = []
        self.best_groups = []               # a list of best contact lists
        self.best_group_scores = []                 # a list of best contact lists
        
        self.FANCY_BOUNDS = fancy_bounds    # flag for bound condition
                                           
        self.RESTONLY = restonly            # flag for only considering incompatibility
                                            # with restraints (not other contacts)
        
        self.counter = 0
        self.WRITE_COUNTER_EVERY = 2500 
        
        self.VERBOSE = False
        self.MAX_BEST_GROUPS = 100      
        self.MAX_ITERS = 100000

        ### Put contacts in order from largest to smallest punch
        for c,p in self.pdict.iteritems():
            if len(self.contacts)>0:
                i=0
                while (i<len(self.contacts)):
                    if (p > self.weights[i]):
                        self.contacts.insert(i,c)
                        self.weights.insert(i,p)
                        break
                    i=i+1
            else:
                 self.contacts.append(c)
                 self.weights.append(p)
            
        # initialize the selection list with 0....0
        for i in range(0,len(self.contacts)):
             self.select.append(0)
        
        # build up the incompatibility dict
        for i in range(0, len(self.contacts) ):
            self.incompatible[ self.contacts[i] ] = []
        for i in range(0, len(self.contacts) ):
          ci = self.contacts[i]
          if self.RESTONLY == False:
            for j in range(i+1, len(self.contacts) ):
              cj = self.contacts[j]
            
              if self.too_close( ci, cj ):
                self.incompatible[ci].append( cj )
                self.incompatible[cj].append( ci )

          for j in range(0,len(self.restrained_contacts)):
            cj = self.restrained_contacts[j] 
            if self.too_close( ci, cj ):
                self.incompatible[ci].append( cj )
                             

    def __repr__( self ):
    
        tmp = 'contact\tweight\tincompatible\n'
        for i in range(0,len(self.contacts)):       
            c = self.contacts[i]
            tmp = tmp + repr( c )+( '\t%1.4f'%self.weights[i] )+'\t'+repr( self.incompatible[c] )+'\n'
        
        return tmp


    def too_close( self, ci, cj ):        
        if ( abs(ci[0]-cj[0]) + abs(ci[1]-cj[1]) ) < self.incompatible_cutoff:
            return True
        else:
            return False
           
    
    def explore( self, nlimit, selectmax):

        self.counter = self.counter + 1
        if self.counter > self.MAX_ITERS:
            return


        if self.VERBOSE:
          if (self.counter % self.WRITE_COUNTER_EVERY ) == 0:
            print self.counter, 'cgroup.explore() iters' 
        
        if selectmax == len(self.contacts):
           return
            
        self.select[selectmax] = 1

        if DEBUG:
          print self.select
        
        # First, see if this new selected contact is kosher
        newcontact = self.contacts[selectmax]
        
        # 1) Is the newly selected contact compatible with 
        #    chosen contacts so far and restrained contacts?
        compatible = 1
        if self.RESTONLY == False:
          for i in range(0,selectmax):
            ci = self.contacts[i]
            if self.select[i] == 1:
              if self.incompatible[newcontact].count( ci ) > 0:
                  compatible = 0
                  break 
        for i in range(0,len(self.restrained_contacts)):
          ci = self.restrained_contacts[i]
          if self.incompatible[newcontact].count( ci ) > 0:
              compatible = 0
              break     
                        
        if compatible == 0:
            if DEBUG:
              print 'contact',self.contacts[i],'is incompatible'
              print
            
           
        # 2)  Have we reached the bound condition?
        score = 0.0
        possiblescore = 0.0
        for i in range(0,selectmax+1):
            score = score + self.select[i]*self.weights[i]
        for i in range(selectmax+1, len(self.weights)):
            accept = 1
            if self.FANCY_BOUNDS:       
              for j in range(0,selectmax+1):
                if self.incompatible[self.contacts[i]].count( self.contacts[j] ) > 0:
                  accept = 0
                  break
                                
            if accept:
                possiblescore = possiblescore + self.weights[i]
            
        if (score+possiblescore) < self.best_score:
            if DEBUG:
              print '*over bound!*'
              print 'score\tpossiblescore\tbestscore'
            self.select[selectmax] = 0
            return 'done'
            
        # 3) Have we reached the top N limit in number of contacts to be considered?
        if sum(self.select[0:selectmax+1]) > nlimit:
            if DEBUG:
              print 'Over the contact limit!'
              print
            self.select[selectmax] = 0
            return 'done'


        if DEBUG:
          print 'selectmax', selectmax, 'score:',score, 'best:',self.best_score
          print 
        
        if compatible: 
          if score > self.best_score:
            self.best_score = score
            self.best_select = []
            self.best_contacts =[]
            for i in range(0,len(self.select)):
                self.best_select.append( self.select[i] )
                if self.select[i] == 1:         
                    self.best_contacts.append( self.contacts[i] )
            if self.VERBOSE:
              print 'New best score:', self.best_score
              print '\t',self.best_contacts

          ilastbest = min( len(self.best_group_scores), self.MAX_BEST_GROUPS ) - 1
          if ilastbest >= 0:
           if score > self.best_group_scores[ilastbest]:
            for i in range(0,ilastbest+1):
              if score > self.best_group_scores[i]:
                self.best_groups.insert(i, [] )
                for j in range(0,len(self.select)):
                  if self.select[j] == 1:               
                    self.best_groups[i].append( self.contacts[j] )
                if len(self.best_groups) > self.MAX_BEST_GROUPS:
                  self.best_groups.pop()
              
                self.best_group_scores.insert(i, score )
                if len(self.best_group_scores) > self.MAX_BEST_GROUPS:
                  self.best_group_scores.pop()
                break
          else:
                self.best_groups.insert(0, [] )
                for j in range(0,len(self.select)):
                  if self.select[j] == 1:               
                      self.best_groups[0].append( self.contacts[j] )
                self.best_group_scores.insert(0, score )

          ## If weve gotten this far, then try to add another contact
          self.explore( nlimit, selectmax+1 )
                

        ## After we;ve reached the end of the explore() branch, try
        ## skipping this contact (0), and adding another next-in-line one (1)
        self.select[selectmax] = 0
        self.explore( nlimit, selectmax+1 ) 
        
        return




########################
### Functions

def bestn( pdict, n, eco, restrained_contacts, fancy_bounds, restonly ):

    ### return a ranked list of the set N contacts with the greatest score sum

    cgroup = ContactGroup( pdict, eco, restrained_contacts, fancy_bounds, restonly)
    if fancy_bounds == True:
        print '\nFinding the best', n, 'contacts by: BRANCH-AND-BOUND...'
    else:
        print '\nFinding the best', n, 'contacts by: ENUMERATION...'
    cgroup.explore( n, 0)

    return [ cgroup.best_score, cgroup.best_contacts, cgroup.best_groups, cgroup.best_group_scores ]



