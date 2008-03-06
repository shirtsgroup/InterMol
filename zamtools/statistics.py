#!/usr/bin/env python

#LAST MODIFIED: 02-11-07


from numpy import *
import scipy, scipy.optimize
import copy


def LogisticRegression(x, y, Verbose = True, cInit = None):
  """Returns an array of coefficients, the overall logprob,
and the logprobs for each data point..
x[N,M] are the explanatory variables and y[N] are the binary
variables describing the classification."""
  N, M = x.shape
  if cInit is None:
    c = zeros(M+1, float)
  else:
    c = copy.deepcopy(cInit)
  #make a probability function for a single obs
  def logprob(xi, c):
    t = c[0] + dot(xi, c[1:])
    if t < 0:
      return t - log(1. + exp(t))
    else:
      return -log(1. + exp(-t)) 
  #make a liklihood function for a single obs
  def fobs(xi, yi, c):
    t = c[0] + dot(xi, c[1:])
    if t < 0:
      s = -log(1. + exp(t))
      return yi * (t + s) + (1. - yi) * s
    else:
      s = -log(1. + exp(-t))
      return yi * s + (1. - yi) * (s - t)      
  #make an overall likelihood function
  def fobj(c):
    return -sum([fobs(x[i], y[i], c) for i in range(N)])
  #optimize
  c = scipy.optimize.fmin(fobj, c, disp = int(Verbose))
  return c, -fobj(c), array([logprob(x[i], c) for i in range(N)], float)

def FindExplanatory(x, y, Labels = None):
  """Performs an analysis to determine the most explanatory variables."""
  N, M = x.shape
  if Labels is None: Labels = ["Var%d" % i for i in range(M)]
  VarList1 = []
  Results = []
  VarList2 = range(M)
  cInit = [0., 0.]
  while len(VarList2) > 0:
    Scores = []
    for VarInd in VarList2:
      NewVarList = VarList1 + [VarInd]
      Newx = x.take(NewVarList, axis=1)
      Coefs, LogProb, LogProbList = LogisticRegression(Newx, y, Verbose = False,
                                                       cInit = array(cInit, float))
      Scores.append((LogProb, VarInd, Coefs))
      print NewVarList, LogProb
    Scores.sort()
    LogProb, VarInd, Coefs = Scores[-1]
    cInit = list(Coefs) + [0.]
    VarList1.append(VarInd)
    Results.append((copy.deepcopy(VarList1), Coefs, LogProb))
    VarList2 = [i for i in VarList2 if not i == VarInd]
    ThisLabels = ", ".join([Labels[i] for i in VarList1])
    ThisCoefs = ", ".join(["%.4f" % c for c in Coefs])
    print ThisLabels, '\n', ThisCoefs, '\n', LogProb, '\n'


  
  
    
