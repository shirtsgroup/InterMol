#!/usr/bin/env python

#LAST MODIFIED: 04-10-07

from numpy import *
import os, sys, copy
import sequence


#GLOBALS
MinPropRest = 1
MaxPropRest = 3

#these are coefficients for (fraglength, maxCO, prior, cprob, pmfscore, punchstab, punchcoop, meso);
#fraglength must go from low to high; maxCO must go from low to high (or 0 for all);
#next is prior probabilities [(maxco, prior), ...]
LogitParams = [(8,  0, -2.5916,  2.9295,  0.0000, -0.0820,  0.0000, -0.0787),
               (12, 0, -2.3686,  2.5255,  0.0000, -0.0373,  0.0000,  0.1154),
               (16, 0, -1.9864,  2.7035,  0.0107, -0.0262,  0.0000,  0.0306)]

#These are with contact order (doesn't work so well):
##LogitParams = [(8,  4, -0.8016,  0.5690,  0.0000,  0.0000,  0.0000,  0.0000),
##               (8,  0, -3.6160,  2.2618,  0.0000, -0.0322,  0.0128,  0.0000),
##               (12, 4, -0.6491,  0.2837,  0.0602,  0.0000,  0.0000,  0.0000),
##               (12, 0, -3.0935,  2.0215,  0.0000, -0.0334,  0.0000,  0.0000),
##               (16, 4,  0.5020,  1.7225,  0.0000, -0.0275,  0.0000, -0.3364),
##               (16, 0, -3.0505,  0.6961,  0.0000,  0.0000,  0.0000,  0.0000)]


#these add priors and correct for correlations among different measurements;
#these are coefficients for (maxCO, b0, b1, b2, b3, b4)
CorrelParams = [(4, 0.9171, 0.0337, 0.9892, 0.0000, 1.0000, 0.0000),
                (0, 0.9171, 0.0337, 0.9892, 0.0000, 1.0000, 0.0000)]





#======== PRIOR PROBABILITY FUNCTIOBS ========

def AvgNNatCont(NRes):
  "Returns the average number of CO >= 3 native contacts for given length."
  return NRes * (4.2607 * NRes - 2.6029) / (NRes + 41.33)

def AvgNTotCont(NRes):
  "Returns the average total number of CO >= 3 contacts for given length."
  return NRes*(NRes - 5.) / 2. + 3.

def AvgLogit(NRes):
  "Returns the average ratio for PNat/PNon."
  x = AvgNNatCont(NRes) / AvgNTotCont(NRes)
  return x / (1. - x)

def AvgLogPNat(NRes):
  "Returns the average native probability for a given length."
  return log(AvgNNatCont(NRes) / AvgNTotCont(NRes))

def AvgLogPCont(CO):
  "Returns the average contact probability as a function of contact order."
  return -3.52248 + 584.59511 * (CO - 3.78480) / (CO**4.10343 - 249.61419)


#======== PROBABILITY MANIPULATIORS ========

def CalcLogit(CProb, PMFScore, PunchStab, PunchCoop, MesoEntropy, NRes, CO):
  "Returns the log ratio of P(Nat)/P(Non) for a contact probability."
  for (n, c, b0, b1, b2, b3, b4, b5) in LogitParams:
    if n >= NRes and (c >= CO or c == 0): break
  return b0 + b1*CProb + b2*PMFScore + b3*PunchStab + b4*PunchCoop \
          + b5*MesoEntropy + log(AvgLogit(NRes) / AvgLogit(n))

def LogitToLogP(Logit):
  "Converts a log ratio to a log probability."
  x = where(Logit < 0, Logit, 0.)
  y = where(Logit > 0, -Logit, 0.)
  return x - log(exp(x) + exp(y))
  
def LogPToLogit(LogP):
  "Converts a log prob to a log ratio."
  return LogP - log(1. - exp(LogP))

def LogitToP(Logit):
  "Converts a log ratio to a probability."
  return exp(LogitToLogP(Logit))

def CalcScore(LogitSum, Count, CO):
  "Calculates scores and removes correlated measurements."
  for (c, b0, b1, b2, b3, b4, b5) in CorrelParams:
    if c >= CO or c == 0: break
  if Count > 0:
    Ratio = (Count - 1.) / Count
  else:
    Ratio = 0
  return LogitSum * (1. - b2 * Ratio) + b0 + b1*Count + b3 * CO / (CO + b4)


#======== RESTRAINT CHOOSING ========

def ChooseRest(Scores, MinRest = 1, MaxRest = 3, MinSimilarity = 4.5, MinCO = 3,
               SimilarityCOCoef = 0):
  """Returns a list of restraints for a scoring matrix."""
  l = []
  N = len(Scores)
  Scores = Scores.copy()
  DScore = Scores.max() - Scores.min()  
  #black out minimum CO
  for i in range(N):
    for j in range(max(0, i-MinCO+1), min(N, i+MinCO)):
      Scores[i,j] = 1.e-200
  #start picking restraints
  while len(l) < MaxRest and len(l) < N:
    ind = argmax(Scores)
    a, b = int(ind / N), int(ind % N)
    l.append((Scores[a,b], a, b))
    #black out the scores and any surrounding scores
    for i in range(N):
      for j in range(i, N):
        #similarity is defined in terms of distance and CO
        CODiff = abs(abs(a-b) - abs(i-j))
        Similarity = abs(a - i) + abs(b - j) + SimilarityCOCoef * CODiff
        if Similarity < MinSimilarity:
          Scores[i,j] -= DScore
          Scores[j,i] -= DScore
  return l  


#========CLASSES========

class ContactDatabaseClass:

  def __init__(self, zd = None, Update = True):
    """Initializes the class using a particular protein sequence."""
    if zd is None:
      self.N = 0
      self.Seq = []
    else:
      self.N = len(zd.Seq)
      self.Seq = copy.deepcopy(zd.Seq)
    self.Clear()
    #update the database
    if not zd is None and Update:
      for Frag in zd.GetDoneFrags(): self.UpdateFrag(zd, Frag)
      self.UpdateScores()

  def Clear(self):
    """Clears all stored data."""
    self.LogitSum = zeros((self.N,self.N), float)
    self.Scores = zeros((self.N,self.N), float)
    self.Counts = zeros((self.N,self.N), float)
    CO = fromfunction(lambda i,j: abs(i-j), (self.N,self.N))

  def UpdateContact(self, a, b, Logit):
    """Updates contact probability ratios based on log P(Nat)."""
    self.LogitSum[a,b] += Logit
    self.LogitSum[b,a] += Logit
    self.Counts[a,b] += 1.
    self.Counts[b,a] += 1.

  def UpdateScores(self):
    """Updates all scores from contact probability data."""
    for i in range(self.N):
      for j in range(i+1, self.N):
        CO = j - i
        self.Scores[i,j] = CalcScore(self.LogitSum[i,j], self.Counts[i,j], CO)
        self.Scores[j,i] = self.Scores[i,j]

  def UpdateFrag(self, zd, Frag):
    """Updates contact data, proposed restraints, and fragment score."""
    Frag.LoadAnalData()
    #check for no data
    if len(Frag.CutProb) == 0 or len(Frag.PmfScore) == 0 \
       or len(Frag.PunchCoop) == 0 or len(Frag.PunchStab) == 0:
      Frag.ClearAnalData()
      return
    #loop through contacts
    NRes = Frag.StopRes - Frag.StartRes + 1
    Frag.FragScore = 0.
    Scores = ones((NRes,NRes), float) * 1.e-200
    for (i, (a,b)) in enumerate(Frag.PairList):
      Logit = CalcLogit(Frag.CutProb[i], Frag.PmfScore[i],
                        Frag.PunchCoop[i], Frag.PunchStab[i],
                        Frag.MesoEntropy, NRes, abs(a-b))
      self.UpdateContact(a, b, Logit)
      Frag.FragScore += LogitToP(Logit)
      if not sequence.SaltPair(self.Seq[a], self.Seq[b]) and not (a,b) in Frag.PriorRest:
        Scores[a - Frag.StartRes, b - Frag.StartRes] = Logit
        Scores[b - Frag.StartRes, a - Frag.StartRes] = Logit
    Frag.PropRest = ChooseRest(Scores, MinRest = MinPropRest, MaxRest = MaxPropRest)
    Frag.PropRest = [(s, a+Frag.StartRes, b+Frag.StartRes) for (s,a,b) in Frag.PropRest]
    Frag.ClearAnalData()

  def UpdateFrags(self, zd, FragList):
    """Updates a list of fragments."""
    for Frag in FragList: self.UpdateFrag(zd, Frag)

  def WriteScores(self, Filename):
    """Writes scores to a file."""
    s = [["%10.4e" % x for x in y] for y in self.Scores]
    s = "\n".join([" ".join(y) for y in s])
    file(Filename, "w").write(s)

  def WriteHtml(self, Filename):
    """Writes scores to an html file."""
    MaxS = max(max(self.Scores[self.Counts > 0]),  0.001)
    MinS = min(min(self.Scores[self.Counts > 0]), -0.001)
    Frac = arange(5.) / 4.
    s = """<html>
<head>
  <meta content="text/html; charset=ISO-8859-1" http-equiv="content-type">
  <title></title>
  <STYLE TYPE="text/css">
  <!--
  TD
     {
     font-size:8pt;
     font-family:courier;
     }
  -->
  </STYLE>
</head>
<body>
 <small>

 <table style="text-align: left; width: 400px;" border="0" cellpadding="0" cellspacing="1">
  <tbody>
    <tr>
"""
    s += """      <td>score</td>
"""
    for x in Frac[::-1]:
      s += """      <td>%.2f</td>
""" % (x * abs(MinS),)
    for x in Frac[1:]:
      s += """      <td>%.2f</td>
""" % (x * abs(MaxS),)
    s += """    </tr>
    <tr>
      <td>key</td>
"""
    for x in Frac[::-1]:
      v = int(x * 255.)
      s += """      <td style="background-color: rgb(%d, %d, %d);"></td>
""" % (255, 255 - v, 255 - v)
    for x in Frac[1:]:
      v = int(x * 255.)
      s += """      <td style="background-color: rgb(%d, %d, %d);"></td>
""" % (255 - v, 255 - v, 255)
    s += """    </tr>
  </tbody>
 </table>
<br>
"""
    s += """
 <table style="text-align: left; width: 400px;" border="0" cellpadding="0" cellspacing="1">
  <tbody>
    <tr>
"""
    Seq = sequence.SeqToAA1(self.Seq)
    s += """      <td></td>
"""
    for i in range(self.N):
      s += """      <td>%s%-3d</td>
""" % (Seq[i], i+1)
    for i in range(self.N):
      s += """    </tr>
    <tr>
      <td>%s%-3d</td>
""" % (Seq[i], i+1)
      for j in range(self.N):
        if self.Counts[i,j] == 0:
          s += """      <td style="background-color: rgb(%d, %d, %d);"></td>
""" % (200, 200, 200)          
        elif self.Scores[i,j] > 0:
          v = int(255. * abs(self.Scores[i,j] / MaxS))
          s += """      <td style="background-color: rgb(%d, %d, %d);"></td>
""" % (255 - v, 255 - v, 255)
        else:
          v = int(255. * abs(self.Scores[i,j] / MinS))
          s += """      <td style="background-color: rgb(%d, %d, %d);"></td>
""" % (255, 255 - v, 255 - v)
      s += """<td>%s%-3d</td>
""" % (Seq[i], i+1)
    s += """    </tr>
    <tr>
"""
    s += """      <td></td>
"""
    for i in range(self.N):
      s += """      <td>%s%-3d</td>
""" % (Seq[i], i+1)
    s += """    </tr>
  </tbody>
 </table>
<br>

 </small>
</body>
</html>
"""
    file(Filename, "w").write(s)    

