import sys,os,tempfile, string

class Status(object):

  def __init__(self): 
    """A class to keep track of the status of of the simulation system."""

    # logffile lines
    self.loglines = []     # append (cmd, output)
    
  
        
    # Benchmarks
    self.groprepDone = False
    self.genboxDone = False
    self.genionDone = False
    self.minimizeDone = False
    self.equilibrateDone = False
    self.simulateDone = False
  
  
  def writeHistory(self, logfile):
    """Writes logfile lines to file in a standard file format"""
   
    fout = open(logfile,'w')
    for line in self.loglines:
      fout.write(line[0]+'\n')
    fout.close()

  def writeLog(self, logfile):
    """Writes logfile lines to file in a standard file format"""

    fout = open(logfile,'w')
    for line in self.loglines:
      fout.write(line[1]+'\n')
    fout.close()
    

  def __repr__(self): 
    """A class to keep trackk of the status of of the simulation system."""
    
    return "wowwiw!"
