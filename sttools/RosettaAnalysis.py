
import os
import os.path
import shutil
import sys

from EnergyHistogram import *

class RosettaAnalysis:
  def __init__(self, dirName, numIter, numRun, numStruct, numT, myT=[]):
    self.dirName = dirName         # directory name
    self.numIter = numIter         # num iterations
    self.numRun = numRun           # num runs per iter
    self.numStruct = numStruct     # num structure per run
    self.numT = numT               # num temperatures

    self.funcList = []             # list of functions to be called

    self.scoreFn = "scores.sc"     # score file name
    self.decoyFn = "aa1d3z.out"    # decoy file name

    self.numExploreWholeSpace = 0     # num decoys visited each temp
    self.numNotExploreWholeSpace = 0  # num decoys didn't visit each temp
    self.attemptSwapUp = []        # num attempts to swap from each temp to next higher one
    self.acceptSwapUp = []         # num attempts that made it
    self.attemptSwapDown = []      # num attempts to swap from each temp to next lower one
    self.acceptSwapDown = []       # num attempts made it
    for i in range(0, self.numT):
      self.attemptSwapUp.append(0)
      self.acceptSwapUp.append(0)
      self.attemptSwapDown.append(0)
      self.acceptSwapDown.append(0)

    self.scores = []               # num iter by num run by num struct table of scores
    self.rmsds = []                # num iter by num run by num struct table of rmsd
    for iter in range(0, self.numIter):
      self.scores.append([])
      self.rmsds.append([])
      for run in range(0, self.numRun):
        self.scores[iter].append([])
        self.rmsds[iter].append([])

    # paths to perl scripts for running rosetta
    self.rosetta="../../../../source/rosetta.macdebug"
    self.rosettaAB="../../../../scripts/bin/rosettaAB.pl"
    self.cluster="../../../../scripts/bin/cluster.pl"
    self.extract="../../../../scripts/bin/extract.pl"

    self.mainDir = os.path.abspath("./")   # path to directory ran script in

    self.extraStructs = []         # extra structs output if accept final structure
                                   # elements in form [iter, run, score, rmsd]

    self.myT = myT                 # temp list
    self.myWeights = []            # weights for ST

    self.tempSortScore = []        # scores sorted by temp index for reweighting
    for T in range(0, self.numT):
      self.tempSortScore.append([])

  # set the list of functions to be run on each trajectory
  def setFuncList(self, funcList):
    self.funcList = funcList

  # execute each function in funcList on each trajectory
  def execute(self):
    for iter in range(0, self.numIter):
      for run in range(0, self.numRun):
        curDir = os.path.join(self.dirName, "iter" + str(iter), "Run" + str(run))
        for func in self.funcList:
          func(curDir, iter, run)

  # get score info for const T runs
  def analyzeConstTScores(self, runDir, iter, run):
    # should only be one iter
    if self.numIter != 1:
      print "ERROR: should be a single iteration for const T analysis."
      sys.exit(0)

    # should be one run per temp
    if self.numRun != self.numT:
      print "ERROR: number runs should equal number temps."
      sys.exit(0)

    fn = os.path.join(runDir, self.scoreFn)
    f = open(fn, 'r')

    for line in f:
      line = line.strip()
      score = float(line)
      self.scores[iter][run].append(score)
      self.tempSortScore[run].append(score)

    f.close()

  # get weights
  def getWeights(self, bias=0):
    hists = []
    for i in range(0, self.numT):
      hists.append(EnergyHistogram())
      hists[i].set_temp(self.myT[i])
      hists[i].read_values(self.tempSortScore[i])

    # look at each temp pair
    cutoff = 0.001
    finalWeights = [0]
    for i in range(0, self.numT-1):
      Pij = hists[i].acceptRatio(hists[i+1])
      Pji = hists[i+1].acceptRatio(hists[i])

      # how much greater want Pji than Pij
      dP = Pji - Pij - bias

      # use Newton's method to get weights
      while abs(dP) > cutoff:
        dP_prime = -Pji-Pij
        hists[i+1].g = hists[i+1].g - dP/dP_prime

        Pij = hists[i].acceptRatio(hists[i+1])
        Pji = hists[i+1].acceptRatio(hists[i])
        dP = Pji - Pij - bias

      finalWeights.append(hists[i+1].g)
      print "Pij " + str(i) + "->" + str(i+1) + " " + str(Pij)
      print "  Pji " + str(Pji)

    self.myWeights = finalWeights
    print "Final weights " + str(finalWeights)

  # print ST params for runs starting from each temp
  def printSTParams(self):
    for startInd in range(0, self.numT):
      fn = "st_T" + str(startInd) + ".in"
      f = open(fn, 'w')

      f.write("nT: " + str(self.numT) + "\n")
      f.write("ind: " + str(startInd) + "\n")
      f.write("seed: -1" + "\n")
      f.write("nst: 50" + "\n")
      TLine = "T: "
      WLine = "W: "
      for i in range(0, self.numT-1):
        TLine += str(self.myT[i]) + " "
        WLine += str(self.myWeights[i]) + " "
      TLine += str(self.myT[self.numT-1]) + "\n"
      WLine += str(self.myWeights[self.numT-1]) + "\n"
      f.write(TLine)
      f.write(WLine)

      f.close()

  # analyzes score files, pulling out all temp info
  def analyzeScoreFile(self, runDir, iter, run):
    seenTMin = False
    seenTMax = False

    fn = os.path.join(runDir, self.scoreFn)
    f = open(fn, 'r')

    line = f.readline()
    while line:
      line = line.split()
      if "Step:" in line:
        step = int(line[1])
        line = f.readline().split()
        curInd = int(line[1])
        line = f.readline().split()
        curT = float(line[1])
        line = f.readline().split()
        score = float(line[1])
        line = f.readline().split()
        newInd = int(line[1])
        line = f.readline().split()
        newT = float(line[1])
        line = f.readline().split()
        dswap = float(line[1])
        line = f.readline().split()
        pswap = float(line[1])
        line = f.readline().split()
        coinflip = float(line[1])
        line = f.readline().split()
        swap = line[1]

        if newInd<curInd:
          self.attemptSwapDown[curInd] += 1
          if swap=="true":
            self.acceptSwapDown[curInd] += 1
        if curInd<newInd:
          self.attemptSwapUp[curInd] += 1
          if swap=="true":
            self.acceptSwapUp[curInd] += 1

        if curInd == 0:
          seenTMin = True
        if curInd == self.numT-1:
          seenTMax = True

      line = f.readline()

    f.close()

    if seenTMin and seenTMax:
      self.numExploreWholeSpace += 1
    else:
      self.numNotExploreWholeSpace += 1

  # prints temperature information
  def printTempInfo(self):
    print "Num attempts swap up: " + str(self.attemptSwapUp)
    print "Num accepted swap up: " + str(self.acceptSwapUp)
    print "Num attempts swap down: " + str(self.attemptSwapDown)
    print "Num accepted swap down: " + str(self.acceptSwapDown)
    print "Num see whole temp space: " + str(self.numExploreWholeSpace)
    print "Num not see whole temp space: " + str(self.numNotExploreWholeSpace)

  # analyzes trajectories getting scores and rmsds
  def analyzeDecoys(self, runDir, iter, run):
    fn = os.path.join(runDir, self.decoyFn)
    f = open(fn, 'r')

    # skip first two lines (header)
    f.readline()
    f.readline()

    numStructFound = 0
    for line in f:
      line = line.split()

      # only care about score lines for decoys I output
      # ignore final decoy if is one (will start with "S_" since will have been accepted 
      # final struct, all others start with "F_" since I flagged them as unaccepted)
      if "SCORE:" in line[0] and not "S_" in line[len(line)-1]:
        numStructFound += 1
        score = float(line[1])
        rms = float(line[16])

        self.scores[iter][run].append(score)
        self.rmsds[iter][run].append(rms)

      # store info for accepted structures (start with "S_")
      if "SCORE:" in line[0] and "S_" in line[len(line)-1]:
        score = float(line[1])
        rms = float(line[16])

        self.extraStructs.append([iter, run, score, rms])

    f.close()

    if numStructFound != self.numStruct:
      print "ERROR: found " + str(numStructFound) + " decoys in " + fn

  # print rmsd and scores (in two columns) for plotting
  def printRMSDVsScore(self):
    fn = self.dirName + "_RMSD_vs_score.dat"
    f = open(fn, 'w')

    for iter in range(0, self.numIter):
      for run in range(0, self.numRun):
        for struct in range(0, self.numStruct):
          rms = self.rmsds[iter][run][struct]
          score = self.scores[iter][run][struct]

          line = str(rms) + " " + str(score) + "\n"
          f.write(line)

    f.close()

  # like printRMSDVsScore but only prints info for decoy with best score/rmsd
  # best score alway first decoy since sorted by score
  # have to search for best rmsd
  def printBestRMSDVsScore(self):
    fn = self.dirName + "_best_RMSD_vs_score.dat"
    f = open(fn, 'w')

    for iter in range(0, self.numIter):
      for run in range(0, self.numRun):
        # info for best score
        struct = 0
        rms = self.rmsds[iter][run][struct]
        score = self.scores[iter][run][struct]
        line = str(rms) + " " + str(score) + "\n"
        f.write(line)

        # info for best rmsd
        minRMSD = 9999
        minRMSDInd = -1
        for struct in range(0, self.numStruct):
          rms = self.rmsds[iter][run][struct]
          score = self.scores[iter][run][struct]
          if rms < minRMSD:
            minRMSD = rms
            minRMSDInd = struct
        rms = self.rmsds[iter][run][minRMSDInd]
        score = self.scores[iter][run][minRMSDInd]
        line = str(rms) + " " + str(score) + "\n"
        f.write(line)

    f.close()

  # print info on extra accepted structs
  def printExtraAcceptedStructs(self):
    fn = self.dirName + "_extra_RMSD_vs_score.dat"
    fnList = self.dirName + "_extra_list_RMSD_vs_score.dat"
    f = open(fn, 'w')
    f2 = open(fnList, 'w')

    for struct in range(0, len(self.extraStructs)):
      rms = self.extraStructs[struct][3]
      score = self.extraStructs[struct][2]

      line = str(rms) + " " + str(score) + "\n"
      f.write(line)

      iter = self.extraStructs[struct][0]
      run = self.extraStructs[struct][1]
      line = str(iter) + " " + str(run) + " " + line
      f2.write(line)

    f.close()
    f2.close()


if __name__ == "__main__":
  numRuns = 5
  numIter = 1

  # analyze STRuns
  myRA = RosettaAnalysis("STRuns", numIter, numRuns, 1440, 6)
  funcList = [myRA.analyzeScoreFile, myRA.analyzeDecoys]
  myRA.setFuncList(funcList)
  myRA.execute()
  myRA.printTempInfo()
  myRA.printRMSDVsScore()
  myRA.printBestRMSDVsScore()
  myRA.printExtraAcceptedStructs()

  # analyze nonSTRuns
  myRA = RosettaAnalysis("nonSTRuns", numIter, numRuns, 1440, 6)
  funcList = [myRA.analyzeDecoys]
  myRA.setFuncList(funcList)
  myRA.execute()
  myRA.printRMSDVsScore()
  myRA.printBestRMSDVsScore()
  myRA.printExtraAcceptedStructs()


