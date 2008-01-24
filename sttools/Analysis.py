
import os
import os.path
import shutil
import sys

from EnergyHistogram import *

numT = 56
pNum = "3651"
myT = [285, 288, 291, 294, 297, 300, 303, 306, 309, 313, 317, 321, 325, 329, 333, 337, 341, 345, 349, 353, 357, 361, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 426, 432, 438, 444, 450, 456, 462, 468, 475, 482, 489, 496, 503, 510, 517, 525, 533, 541, 549, 557, 566, 575]

allEnergy = []
hists = []

for i in range(0, numT):
  print "Read T " + str(myT[i]) + "..."
  allEnergy.append([])
  fn = "p" + pNum + ".T." + str(i) + ".P"
  f = open(fn, 'r')
  for line in f:
    val = float(line.strip())
    allEnergy[i].append(val)
  f.close()

  hists.append(EnergyHistogram())
  hists[i].set_temp(myT[i])
  hists[i].read_values(allEnergy[i])

cutoff = 0.001
bias = 0
finalWeights = [0]
for i in range(0, numT-1):
  Pij = hists[i].acceptRatio(hists[i+1])
  Pji = hists[i+1].acceptRatio(hists[i])

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

print "Final weights " + str(finalWeights)


