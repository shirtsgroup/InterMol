#!/usr/bin/env python


#DESC: Contains non-class functions for protein.py


#LAST MODIFIED: 03-11-07

import sys

from numpy import *
from geometry import *
from proteinconst import *
import sequence



#======== FORCE FIELD PARAMS ========

def GetFFData(Element, AtomName, ResName):
    """Returns force field parameters (Sigma, SqrtEps, PCharge)"""
    #check aliases
    if not ResName in FFResData:
        aa = sequence.AAInstAlias(ResName)
        if not aa is None: ResName = aa.Name
    #first check residue params
    if ResName in FFResData:
        if AtomName in FFResData[ResName]:
            return FFResData[ResName][AtomName]
    #check atom name
    if AtomName in FFAtomData:
        return FFAtomData[AtomName]
    #check element
    if Element in FFElementData:
        return FFElementData[Element]
    else:
        return FFElementData["*"]


#======== BONDS ========

def GetBonds(Bonds):
    """Returns lists of CO<=2 (1-2 and 1-3) and CO=3 (1-4) bonds."""
    b = sorted([(min(x), max(x)) for x in Bonds])
    Bonds13, Bonds14 = [], []
    for (a1,b1) in b:
        #check a links
        clist = [b2 for (a2,b2) in b if a1==a2 and b2 < b1 and not b1==b2] +\
              [a2 for (a2,b2) in b if a1==b2 and a2 < b1 and not b1==a2]
        Bonds13.extend([(min(c,b1), max(c,b1)) for c in clist])
        clist = [b2 for (a2,b2) in b if a1==a2 and not b1==b2] +\
              [a2 for (a2,b2) in b if a1==b2 and not b1==a2]
        dlist = [a2 for (a2,b2) in b if b1==b2 and not a1==a2] +\
              [b2 for (a2,b2) in b if b1==a2 and not a1==b2]
        Bonds14.extend([(min(c,d), max(c,d)) for c in clist for d in dlist])
    Bonds1213 = b + Bonds13
    Bonds1213.sort()
    Bonds14.sort()
    Bonds1213 = array(Bonds1213, int)
    Bonds14 = array(Bonds14, int)
    return Bonds1213, Bonds14 

def ConvertBonds(Bonds, NAtom):
    """Returns a sortedlist of IDs for bonds, where ID = i*N - (i+1)(i+2)/2 + j
for i < j."""
    ID = [i*NAtom - (i+1)*(i+2)/2 + j for (i,j) in Bonds]
    ID.sort()
    #get unique values
    ID = [x for (i,x) in enumerate(ID) if i == 0 or x != ID[i-1]]
    return array(ID, int)

def GetBondMasks(i, N, Bonds23, Bonds4):
    """Returns an array of indices for bonds."""
    j = i + 1
    MinID = i*N - (i+1)*(i+2)/2 + j
    MaxID = MinID + N - j - 1
    b23 = Bonds23[logical_and(Bonds23 >= MinID, Bonds23 <= MaxID)] - MinID + j
    b4 = Bonds4[logical_and(Bonds4 >= MinID, Bonds4 <= MaxID)] - MinID + j
    MaskNot23 = ones(N, bool)
    MaskNot23[b23] = False
    Mask4 = zeros(N, bool)
    Mask4[b4] = True
    return MaskNot23, Mask4


#======== CONTACT FUNCTIONS ========

def GetContactList(Pos1, Ind1 = None, Pos2 = None, Ind2 = None, 
                   Radius = 8., MinCO = 3):
    """Returns a list of contacts for one or two position matrices.
  Ind1, Ind2 are list of indices for positions in Pos1, Pos2.
  Pos1, Pos2 are arrays of positions.  Radius is the contact radius.
  MinCO is minimum contact order, only used when Pos2 is not present."""
    l = []
    N1 = len(Pos1)
    if Ind1 is None: Ind1 = range(0,N1)
    if not Pos2 is None:
        N2 = len(Pos2)
        if Ind2 is None: Ind2 = range(0,N2)
        for i in range(0, N1):
            Vecs = Pos2 - Pos1[i,:]
            Vecs = (Vecs*Vecs).sum(axis=1)
            ladd = [(Ind1[i], Ind2[x])
                    for x in flatnonzero(Vecs < Radius*Radius)]
            l.extend(ladd)
    else:
        for i in range(0, N1-1):
            Vecs = Pos1[i+1:] - Pos1[i,:]
            Vecs = (Vecs*Vecs).sum(axis=1)
            ladd = [(Ind1[i], Ind1[i+1+x])
                    for x in flatnonzero(Vecs < Radius*Radius)]
            l.extend(ladd)
        #prune for CO
        l = [(a,b) for (a,b) in l if b - a >= MinCO]
    return l

def GetContactMap(Pos1, Ind1 = None, NInd1 = None,
                  Pos2 = None, Ind2 = None, NInd2 = None,
                  Radius = 8., MinCO = 3):
    """Returns a contact map for one or two position matrices.
  Ind1, Ind2 are list of indices for positions in Pos1, Pos2.
  Pos1, Pos2 are arrays of positions.  Radius is the contact radius.
  MinCO is minimum contact order, only used when Pos2 is not present."""
    N1 = len(Pos1)
    if Ind1 is None: Ind1, NInd1 = range(0,N1), N1
    if not Pos2 is None:
        N2 = len(Pos2)
        if Ind2 is None: Ind2, NInd2 = range(0,N2), N2
        Map = zeros((NInd1, NInd2), float)
        for i in range(0, N1):
            Vecs = Pos2 - Pos1[i,:]
            Vecs = (Vecs*Vecs).sum(axis=1)
            Vals = (Vecs < Radius*Radius).astype(float)
            Map[(Ind1[i], Ind2)] = Vals
    else:
        Map = zeros((NInd1, NInd1), float)
        for i in range(0, N1-1):
            Vecs = Pos1[i+1:,:] - Pos1[i,:]
            Vecs = (Vecs*Vecs).sum(axis=1)
            Vals = (Vecs < Radius*Radius).astype(float)
            Map[(Ind1[i], Ind1[i+1:])] = Vals
        #prune for CO
        for i in range(0, N1):
            Map[(i, range(i, min(i+MinCO, N1)))] = 0.
        #update other half of contact map
        Map = Map + Map.transpose()
    return Map


#======== INDEX FUNCTIONS ========

def GetInd(Ind, ResInd = None):
    """Returns new Ind where > 0."""
    Mask = Ind >= 0
    if not ResInd is None:
        ResMask = zeros_like(Mask)
        ResMask.put(ResInd, 1)
        Mask = logical_and(ResMask, Mask) 
    return Ind[Mask], where(Mask)[0]

def GetCommonInd(Ind1, Ind2, ResInd = None):
    """Returns new Ind1 and Ind2 where both are > 0."""
    Mask = logical_and(Ind1 >= 0, Ind2 >= 0)
    if not ResInd is None:
        ResMask = zeros_like(Mask)
        ResMask.put(ResInd, 1)
        Mask = logical_and(ResMask, Mask) 
    return Ind1[Mask], Ind2[Mask], where(Mask)[0]


#======== INFORMATIONAL FUNCTIONS ========

def IsHBond(PosN, PosO, PosC, PosH = None):
    """Returns True if positions meet requirements for a backbone hydrogen bond."""
    if PosH is None:
        if Length(PosN - PosO) < HBondNODist:
            Vec, Ang = GetVecMapping(PosN - PosO, PosO - PosC)
            return (abs(Ang) <= HBondNOCAng)
        else:
            return False
    elif HBondDSSP:
        dON = Length(PosO - PosN)
        dCH = Length(PosC - PosH)
        dOH = Length(PosO - PosH)
        dCN = Length(PosC - PosN)
        if dOH > dON or dOH > dCH or dOH > dCN:
            return False
        else:
            return (1/dON + 1/dCH - 1/dOH - 1/dCN) < HBondDSSPDistCut
    else:
        if Length(PosH - PosO) < HBondHODist:
            Vec, Ang = GetVecMapping(PosN - PosH, PosH - PosO)
            return (abs(Ang) <= HBondNHOAng)
        else:
            return False

def GetRamaProb(Phi, Psi, RamaProb):
    """Returns the probability associated with a Phi, Psi angle."""
    RamaNPhi, RamaNPsi = RamaProb.shape
    RamaDPhi = 360. / RamaNPhi
    RamaDPsi = 360. / RamaNPsi
    PhiInd = int(clip((Phi + 180.) / RamaDPhi, 0, RamaNPhi - 1))
    PsiInd = int(clip((Psi + 180.) / RamaDPsi, 0, RamaNPsi - 1))
    return RamaProb[PhiInd, PsiInd]

def GetRandomDih(RamaProb):
    """Returns a random dihedral pair and its probability."""
    RamaNPhi, RamaNPsi = RamaProb.shape
    RamaDPhi = 360. / RamaNPhi
    RamaDPsi = 360. / RamaNPsi
    r = random.random()
    CumProb = 0.
    for i in range(RamaNPhi):
        for j in range(RamaNPsi):
            CumProb += RamaProb[i,j]
            if r < CumProb:
                Phi = -180. + RamaDPhi * (i + random.random())
                Psi = -180. + RamaDPsi * (j + random.random())
                return Phi, Psi, RamaProb[i,j]



#======== MONTE CARLO MOVES ========

def GetRandomInd(Prob):
    """Returns a random index of 0 to len(Prob)-1, where the
weight of each index is given by Prob[i].  sum(Prob) must equal 1."""
    CumProb = 0.
    r = random.random()
    for i, p in enumerate(Prob):
        CumProb += p
        if r < CumProb: break
    return i


class MCMoveClass:

    def __init__(self, p):
        """Initializes a MC move for a proteinclass object;
To be modified for derived classes."""
        self.InitBase()

    def InitBase(self, D = 1., MinD = 0., MaxD = 1.e200,
                 DExplore = None, DRefine = None, DTweak = None,
                 W = 1., WExplore = 1., WRefine = 1., WTweak = 1.,
                 TargetAcc = 0.5, NMove = 0, DFinal = None):
        """Initializes base parameters describing the MC Move."""
        self.D = D
        #initial and final d
        self.DInit = D
        self.DFinal = DFinal
        if self.DFinal is None: self.DFinal = self.DInit
        #minimum and maximum possible d
        self.MinD = MinD
        self.MaxD = MaxD
        #target acceptance ratio
        self.TargetAcc = TargetAcc
        #number of possible moves for this guy
        self.NMove = NMove
        #typical values of d for mc annealing stages
        self.DExplore = DExplore
        if self.DExplore is None: self.DExplore = MaxD / 4.
        self.DRefine = DRefine
        if self.DRefine is None: self.DRefine = MaxD / 16.
        self.DTweak = DTweak
        if self.DTweak is None: self.DTweak = MaxD / 64.
        #the probability weight of this move
        self.W = W
        #the weights of this move for mc annealing stages
        self.WExplore = WExplore
        self.WRefine = WRefine
        self.WTweak = WTweak
        #reset counters
        self.Reset()

    def SetExplore(self, Weight = True, Reset = True):
        """Sets the max move to an exploration value."""
        self.D = self.DExplore
        if Weight: self.W = self.WExplore
        if Reset: self.Reset()

    def SetRefine(self, Weight = True, Reset = True):
        """Sets the max move to a refining value."""
        self.D = self.DRefine
        if Weight: self.W = self.WRefine
        if Reset: self.Reset()

    def SetTweak(self, Weight = True, Reset = True):
        """Sets the max move to a tweaking/minimizing value."""
        self.D = self.DTweak
        if Weight: self.W = self.WTweak
        if Reset: self.Reset()

    def Reset(self, ResetD = False):
        """Resets counters and D."""
        if ResetD: self.D = self.DInit
        self.NAcc = 0.
        self.NAtt = 0.

    def Progress(self, FracDone = None):
        """Updates the maximum move if it is to be scaled with simulation progress."""
        if not FracDone is None and not self.DInit == self.DFinal:
            self.D = self.DInit * (1. - FracDone) + self.DFinal * FracDone

    def Make(self, p, ScoreFn, OldE, FracDone = None):
        """This to be filled in by derived classes.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False).
To be modified for derived classes."""
        self.NAtt += 1
        if Acc: self.NAcc += 1.       

    def Update(self):
        """Updates the maximum move to achieve target acceptance."""
        if self.NAtt > 0:
            self.D *= (self.NAcc / self.NAtt + 0.1) / (self.TargetAcc + 0.1)
            self.D = min(self.MaxD, max(self.MinD, self.D))

    def Active(self):
        """Indicates whether there are moves available."""
        return self.NMove > 0


class MCPhiPsiMoveClass(MCMoveClass):
    """Makes a random phi or psi perturbation"""

    def __init__(self, p, MaxAng = 20, ResProb = None, CenterChain = 0):
        if ResProb is None: ResProb = ones(len(p), float)
        ResProb = array(ResProb, float)
        ResProbPhi, ResProbPsi = ResProb.copy(), ResProb.copy()
        for i, x in enumerate(p.Seq):
            if x in NoDihedrals: ResProbPhi[i], ResProbPsi[i] = 0., 0.
            if x in FixedPhi: ResProbPhi[i] = 0.
            if x in FixedPsi: ResProbPsi[i] = 0.
        self.ResProb = concatenate((ResProbPhi, ResProbPsi))
        if sum(self.ResProb) > 0:
            self.ResProb = self.ResProb / sum(self.ResProb)
        self.NRes = len(p)
        self.CenterChain = CenterChain
        self.InitBase(D = MaxAng, MinD = 0., MaxD = 180.,
                      DExplore = 45., DRefine = 10., DTweak = 3.6,  #CHANGED
                      WExplore = 1., WRefine = 1., WTweak = 1.,
                      TargetAcc = 0.5, NMove = sum(self.ResProb > 0.))

    def Make(self, p, ScoreFn, OldE, FracDone = None):
        """Makes a random phi/psi move.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)
        #pick a random residue and phi or psi
        self.ResNum = GetRandomInd(self.ResProb)
        self.MoveType = self.ResNum / self.NRes
        self.ResNum = self.ResNum % self.NRes
        self.OldPos = p.Pos.copy()
        dAng = (2.*random.random() - 1.) * self.D
        if self.MoveType == 0:
            p.RotatePhi(self.ResNum, dAng)
        else:
            p.RotatePsi(self.ResNum, dAng)
        if not self.CenterChain is None: p.Center(ChainNum = self.CenterChain)
        NewE = ScoreFn(p, None)
        return NewE, 0.

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False)."""
        self.NAtt += 1
        if Acc:
            self.NAcc += 1.
        else:
            p.Pos = self.OldPos
        self.OldPos = None


class MCRamaMoveClass(MCMoveClass):
    """Makes a hop in the Ramachandran plot according to RamaProb probabilities."""

    def __init__(self, p, RamaProb, ResProb = None, CenterChain = 0):
        if ResProb is None: ResProb = ones(len(p), float)
        ResProb = array(ResProb, float)
        for i, x in enumerate(p.Seq):
            if x in NoDihedrals: ResProb[i] = 0.
        self.ResProb = ResProb
        if sum(self.ResProb) > 0:
            self.ResProb = self.ResProb / sum(self.ResProb)
        self.RamaProb = RamaProb
        self.CenterChain = CenterChain
        self.InitBase(D = 0., MinD = 0., MaxD = 0.,
                      WExplore = 0., WRefine = 0., WTweak = 0.,
                      TargetAcc = 0.2, NMove = sum(ResProb > 0.))

    def Make(self, p, ScoreFn, OldE, FracDone = None):
        """Makes a random phi/psi move.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)
        self.ResNum = GetRandomInd(self.ResProb)
        self.OldPos = p.Pos.copy()
        #get the current angles
        Phi, Psi = p.PhiPsi(self.ResNum)
        OldP = GetRamaProb(Phi, Psi, self.RamaProb)
        #get a random angle
        Phi, Psi, NewP = GetRandomDih(self.RamaProb)
        #rotate and get the final angles and prob
        Phi, Psi = p.RotateToPhiPsi(self.ResNum, Phi, Psi)
        if not self.CenterChain is None: p.Center(ChainNum = self.CenterChain)
        NewE = ScoreFn(p, None)
        NewP = GetRamaProb(Phi, Psi, self.RamaProb)
        return NewE, log(max(OldP, 1.e-100) / max(NewP, 1.e-100))

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False)."""
        self.NAtt += 1
        if Acc:
            self.NAcc += 1.
        else:
            p.Pos = self.OldPos
        self.OldPos = None


class MCChiMoveClass(MCMoveClass):
    """Makes a random perturbation to a chi angle."""

    def __init__(self, p, MaxAng = 20, ResProb = None):
        if ResProb is None: ResProb = ones(len(p), float)
        ResProb = array(ResProb, float)
        #make a list of chi angles and their probabilities
        self.ChiList, self.ChiProb = [], []
        for rn, r in enumerate(p.Res):
            if ResProb[rn] <= 0.: continue
            cl = [(rn,x) for x in range(len(r.ChiAtoms))
                  if not x in FixedChi.get(r.Name, [])]
            self.ChiList.extend(cl)
            self.ChiProb.extend([ResProb[rn]] * len(cl))
        self.ChiProb = array(self.ChiProb, float)
        if len(self.ChiProb) > 0 and sum(self.ChiProb) > 0:
            self.ChiProb = self.ChiProb / sum(self.ChiProb)
        self.InitBase(D = MaxAng, MinD = 0., MaxD = 180.,
                      DExplore = 90., DRefine = 30., DTweak = 10.,
                      WExplore = 0.1, WRefine = 0.25, WTweak = 1.,
                      TargetAcc = 0.5, NMove = sum(self.ChiProb > 0.))

    def Make(self, p, ScoreFn, OldE, FracDone = None):
        """Makes a random phi/psi move.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)
        self.ChiInd = GetRandomInd(self.ChiProb)
        #random chi angle
        dAng = (2.*random.random() - 1.) * self.D
        self.ResNum, self.ChiNum = self.ChiList[self.ChiInd]
        r = p.Res[self.ResNum]
        #save old position and get old energy
        self.a1, self.a2 = r.StartAtom, r.StartAtom + len(r.Atoms)
        self.OldPos = p.Pos[self.a1:self.a2,:].copy()
        dE = -ScoreFn(p, self.ResNum)
        #rotate angle and get new energy
        p.RotateChi(self.ResNum, self.ChiNum, dAng)
        dE += ScoreFn(p, self.ResNum)
        return OldE + dE, 0.

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False)."""
        self.NAtt += 1
        if Acc:
            self.NAcc += 1.
        else:
            p.Pos[self.a1:self.a2,:] = self.OldPos
        self.OldPos = None


class MCChainRotateClass(MCMoveClass):
    """Rotates chains around a common point or their centroid."""

    def __init__(self, p, MaxAng = 20, ChainProb = None, Point = None):
        NChain = len(p.Chains)
        if ChainProb is None: ChainProb = ones(NChain, float)
        ChainProb = array(ChainProb, float)
        self.ChainProb = ChainProb
        if sum(self.ChainProb) > 0:
            self.ChainProb = self.ChainProb / sum(self.ChainProb)
        self.Point = Point
        self.InitBase(D = MaxAng, MinD = 0., MaxD = 180.,
                      DExplore = 20., DRefine = 10., DTweak = 3.6,
                      WExplore = 1., WRefine = 1., WTweak = 1.,
                      TargetAcc = 0.5, NMove = sum(ChainProb > 0.))

    def Make(self, p, ScoreFn, OldE, FracDone = None):
        """Makes a random phi/psi move.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)
        self.ChainNum = GetRandomInd(self.ChainProb)
        #random chi angle
        dAng = (2.*random.random() - 1.) * self.D
        if self.Point is None:
            Point = p.Centroid(ChainNum = self.ChainNum)
        else:
            Point = self.Point
        #save old position and get old energy
        self.a1, self.a2 = p.ChainAtomNums[self.ChainNum:self.ChainNum+2]
        self.OldPos = p.Pos[self.a1:self.a2,:].copy()
        #rotate and get new energy
        p.Rotate(RandVec(3), dAng, Point = Point, ChainNum = self.ChainNum)
        NewE = ScoreFn(p, None)
        return NewE, 0.

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False)."""
        self.NAtt += 1
        if Acc:
            self.NAcc += 1.
        else:
            p.Pos[self.a1:self.a2,:] = self.OldPos
        self.OldPos = None      


class MCChainTranslateClass(MCMoveClass):
    """Translates chains towards or away from a point."""

    def __init__(self, p, MaxDist = 2., ChainProb = None, Point = None):
        NChain = len(p.Chains)
        if ChainProb is None: ChainProb = ones(NChain, float)
        ChainProb = array(ChainProb, float)
        self.ChainProb = ChainProb
        if sum(self.ChainProb) > 0:
            self.ChainProb = self.ChainProb / sum(self.ChainProb)
        self.Point = Point
        if self.Point is None: self.Point = p.Centroid()
        self.InitBase(D = MaxDist, MinD = 0., MaxD = 1.e200,
                      DExplore = 3., DRefine = 1.0, DTweak = 0.3,
                      WExplore = 1., WRefine = 1., WTweak = 1.,
                      TargetAcc = 0.5, NMove = sum(ChainProb > 0.))

    def Make(self, p, ScoreFn, OldE, FracDone):
        """Makes a random phi/psi move.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)
        self.ChainNum = GetRandomInd(self.ChainProb)
        #random distance
        dR = (2.*random.random() - 1.) * self.D
        Vec = p.Centroid(ChainNum = self.ChainNum) - self.Point
        #save old position and get old energy
        self.a1, self.a2 = p.ChainAtomNums[self.ChainNum:self.ChainNum+2]
        self.OldPos = p.Pos[self.a1:self.a2,:].copy()
            #rotate and get new energy
        p.Translate(dR * UnitVec(Vec), ChainNum = self.ChainNum)
        NewE = ScoreFn(p, None)
        #calculate acceptance factor
        OldRSq = dot(Vec, Vec)
        NewRSq = (sqrt(OldRSq) + dR)**2
        return NewE, log(NewRSq / OldRSq)

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False)."""
        self.NAtt += 1
        if Acc:
            self.NAcc += 1.
        else:
            p.Pos[self.a1:self.a2,:] = self.OldPos
        self.OldPos = None


class MCFragInsertClass(MCMoveClass):
    """Makes a random fragment insertion; this does not obey detailed balance."""

    def __init__(self, p, Frags, CenterChain = 0):
        #Frags is a list of (Weight, StartRes, [[Phi0,Psi0], [Phi1,Psi1], ...])
        TotWeight = float(sum([x[0] for x in Frags]))
        self.Frags = []
        CumProb = 0.
        #calculate cumulative probabilities
        ResUsed = zeros(len(p), int)
        for (Weight, StartRes, DihList) in Frags:
            CumProb += Weight / TotWeight
            self.Frags.append((CumProb, StartRes, DihList))
            ResUsed[StartRes:StartRes + len(DihList)] = 1
        #initialize rest
        self.CenterChain = CenterChain
        self.InitBase(NMove = sum(ResUsed))

    def Make(self, p, ScoreFn, OldE, FracDone = None):
        """Makes a random phi/psi move.
Returns the new energy and the log of P(old)/P(new)."""
        self.Progress(FracDone)
        #pick a random fragment
        x = random.random()
        self.FragNum = 0
        #CAN WE SPEED THIS UP?!!
        while self.Frags[self.FragNum][0] < x and self.FragNum < len(self.Frags) - 1:
            self.FragNum += 1
        #make the move
        self.OldPos = p.Pos.copy()
        (Prob, StartRes, DihList) = self.Frags[self.FragNum]
        for (i, (Phi, Psi)) in enumerate(DihList):
            Phi, Psi = p.RotateToPhiPsi(StartRes + i, Phi, Psi)
        if not self.CenterChain is None: p.Center(ChainNum = self.CenterChain)
        NewE = ScoreFn(p, None)
        return NewE, 0.

    def Accept(self, p, Acc):
        """Processes an acceptance (Acc=True) or rejection (Acc=False)."""
        self.NAtt += 1
        if Acc:
            self.NAcc += 1.
        else:
            p.Pos = self.OldPos
        self.OldPos = None
