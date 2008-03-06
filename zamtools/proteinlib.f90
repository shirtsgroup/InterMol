!bond index and iteration functions

integer function GetBondID(i, j, NAtom)
	integer, intent(in) :: i, j, NAtom
	GetBondID = i * NAtom - (i + 1)*(i + 2)/2 + j
end function

subroutine InitBondIter(NextID, Ind, BondID, NAtom, NBondID)
	implicit none
	integer, intent(in) :: NAtom, NBondID
	integer, intent(inout) :: NextID, Ind
	integer, dimension(NBondID), intent(in) :: BondID
	if (NBondID == 0) then
		NextID = NAtom*NAtom
		Ind = 0
	else
		NextID = BondID(1)
		Ind = 1
	endif
end subroutine

logical function IsBondIter(CurID, NextID, Ind, BondID, NAtom, NBondID)
	!this assumes BondID is sorted and this is called from iteration over a loop 
	!where CurID increases over i = 1, NAtom; j = i+1, NAtom
	implicit none
	integer, intent(in) :: NAtom, NBondID
	integer, intent(inout) :: CurID, NextID, Ind
	integer, dimension(NBondID), intent(in) :: BondID
	IsBondIter = .false.
	if (CurID == NextID) then
		if (Ind < NBondID) then
			Ind = Ind + 1
			NextID = BondID(Ind)
		endif
		IsBondIter = .true.
	endif
end function


!======== CENTROID CALCULATIONS ========

subroutine ResPosCentroid(Pos, AtomResNum, NRes, ResPos, NAtom, Dim)
	implicit none
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer, intent(in) :: NRes
	real*8, dimension(NRes,Dim), intent(out) :: ResPos
	integer, dimension(NRes) :: ResNAtom
	integer :: i, r	
	ResPos = 0.
	ResNAtom = 0
	do i = 1, NAtom
		r = AtomResNum(i) + 1
		ResPos(r, :) = ResPos(r, :) + Pos(i,:)
		ResNAtom(r) = ResNAtom(r) + 1
	enddo
	do i = 1, NRes
		ResPos(i,:) = ResPos(i,:) / real(ResNAtom(i))
	enddo
end subroutine


subroutine ResPosCentroidInd(Pos, AtomResNum, NRes, ResInd, ResPos, NAtom, Dim, NResInd)
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer, dimension(NResInd), intent(in) :: ResInd
	real*8, dimension(NResInd,Dim), intent(out) :: ResPos
	integer, dimension(NResInd) :: ResNAtom
	integer, dimension(NRes) :: NewResNum
	integer :: i, r
	ResPos = 0.
	ResNAtom = 0
	NewResNum = 0
	do i = 1, NResInd
		NewResNum(ResInd(i) + 1) = i
	enddo
	do i = 1, NAtom
		r = NewResNum(AtomResNum(i) + 1)
		if (r == 0) cycle
		ResPos(r, :) = ResPos(r, :) + Pos(i,:)
		ResNAtom(r) = ResNAtom(r) + 1
	enddo
	do i = 1, NResInd
		ResPos(i,:) = ResPos(i,:) / real(ResNAtom(i))
	enddo
end subroutine

subroutine CentroidRange(Pos, Atom1, Atom2, Ret, NAtom, Dim)
	implicit none
	integer, intent(in) :: NAtom, Dim, Atom1, Atom2
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(Dim), intent(out) :: Ret
	integer :: i
	Ret = sum(Pos(Atom1+1:Atom2,:), 1) / real(max(Atom2 - Atom1, 1))
end subroutine

subroutine CentroidAll(Pos, Ret, NAtom, Dim)
	implicit none
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(Dim), intent(out) :: Ret
	integer :: i
	Ret = sum(Pos, 1) / real(max(NAtom, 1))
end subroutine

subroutine Center(Pos, CentroidPos, NewPos, NAtom, Dim)
	implicit none
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(Dim), intent(in) :: CentroidPos
	real*8, dimension(NAtom,Dim), intent(out) :: NewPos
	integer :: i
	do i = 1, Dim
		NewPos(:,i) = Pos(:,i) - CentroidPos(i)
	enddo
end subroutine


!======== OVERLAP SCORE ========

real*8 function OverlapScore(Pos, Cutoff, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer :: i, j
	real*8 :: x, CutSq
	real*8, dimension(Dim) :: Posi
	OverlapScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		do j = i+1, NAtom
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			x = 1. / x
			x = x * x * x
			x = x * x
			OverlapScore = OverlapScore + x
		enddo
	enddo
end function

real*8 function OverlapScore2(Pos1, Pos2, Cutoff, NAtom1, NAtom2, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom1, NAtom2, Dim
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	integer :: i, j
	real*8 :: x, CutSq
	real*8, dimension(Dim) :: Posi
	OverlapScore2 = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom1
		Posi = Pos1(i,:)
		do j = 1, NAtom2
			x = sum((Posi - Pos2(j,:))**2)
			if (x > CutSq) cycle
			x = 1. / x
			x = x * x * x
			x = x * x
			OverlapScore2 = OverlapScore2 + x
		enddo
	enddo
end function


!======== RESIDUE CONTACT SCORE ========

real*8 function ResContactScore(Pos, ResID, EpsMat, Sigma, NAtom, Dim, NEpsMat)
	implicit none
	integer, intent(in) :: NAtom, Dim, NEpsMat
	real*8, intent(in) :: Sigma
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NAtom), intent(in) :: ResID
	real*8, dimension(NEpsMat,NEpsMat), intent(in) :: EpsMat
	real*8, parameter :: CutSq = 1.259921049895  !2.**(1./3.)
	real*8, dimension(Dim) :: Posi
	real*8 :: x, E, invSigSq
	integer :: i, j, IDi, IDj
	ResContactScore = 0.
	invSigSq = 1./(Sigma*Sigma)
	do i = 1, NAtom
		Posi = Pos(i,:)
		IDi = ResID(i)
		if (IDi < 0) cycle
		do j = i+1, NAtom
			IDj = ResID(j)
			if (IDj < 0) cycle
			x = invSigSq * sum((Posi - Pos(j,:))**2)
			if (x < CutSq) then
				E = 1.
			else
				x = 1. / x
				x = x * x * x
				E = 4.*(x - x*x)
			endif
			ResContactScore = ResContactScore + EpsMat(IDi+1, IDj+1) * E
		enddo
	enddo
end function



!======== STERIC SCORE ========

real*8 function StericScore(Pos, Radius, Cutoff, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma
	real*8, dimension(Dim) :: Posi
	StericScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		do j = i+1, NAtom
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			x = x * x
			StericScore = StericScore + x
		enddo
	enddo
end function

real*8 function StericScoreRes(Pos, Radius, Cutoff, ResNum, AtomResNum, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma
	real*8, dimension(Dim) :: Posi
	StericScoreRes = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		if (AtomResNum(i) /= ResNum) cycle
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		do j = 1, NAtom
			if (AtomResNum(j)==ResNum .and. j <= i) cycle
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			x = x * x
			StericScoreRes = StericScoreRes + x
		enddo
	enddo
end function

real*8 function StericScoreBond(Pos, Radius, Cutoff, Bond23, Bond4, Scale4, &
	& NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma
	real*8, dimension(Dim) :: Posi
	real*8 :: Scale
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	StericScoreBond = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				Scale = Scale4
			else
				Scale = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			x = x * x
			StericScoreBond = StericScoreBond + Scale*x
		enddo
	enddo
end function

real*8 function StericScoreResBond(Pos, Radius, Cutoff, ResNum, AtomResNum, &
	& Bond23, Bond4, Scale4, NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma
	real*8, dimension(Dim) :: Posi
	real*8 :: Scale
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	logical :: IsResi, IsResj
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	StericScoreResBond = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		IsResi = (AtomResNum(i) == ResNum)
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				Scale = Scale4
			else
				Scale = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			IsResj = (AtomResNum(j) == ResNum)
			if (.not. (IsResi .or. IsResj)) cycle
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			x = x * x
			StericScoreResBond = StericScoreResBond + Scale*x
		enddo
	enddo
end function


real*8 function StericScore2(Pos1, Pos2, Radius1, Radius2, Cutoff, NAtom1, NAtom2, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom1, NAtom2, Dim
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	real*8, dimension(NAtom1), intent(in) :: Radius1
	real*8, dimension(NAtom2), intent(in) :: Radius2
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma
	real*8, dimension(Dim) :: Posi
	StericScore2 = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom1
		Posi = Pos1(i,:)
		Radiusi = Radius1(i)
		do j = 1, NAtom2
			x = sum((Posi - Pos2(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius2(j)
			x = Sigma * Sigma / x
			x = x * x * x
			x = x * x
			StericScore2 = StericScore2 + x
		enddo
	enddo
end function


!======== LENNARD JONES INTERACTIONS ========

real*8 function LJScore(Pos, Radius, SqrtEps, Cutoff, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma, SqrtEpsi
	real*8, dimension(Dim) :: Posi
	LJScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		do j = i+1, NAtom
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			LJScore = LJScore + SqrtEpsi * SqrtEps(j) * x * (x - 2.)
		enddo
	enddo
end function

real*8 function LJScoreRes(Pos, Radius, SqrtEps, Cutoff, ResNum, AtomResNum, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma, SqrtEpsi
	real*8, dimension(Dim) :: Posi
	LJScoreRes = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		if (AtomResNum(i) /= ResNum) cycle
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		do j = 1, NAtom
			if (AtomResNum(j)==ResNum .and. j <= i) cycle
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			LJScoreRes = LJScoreRes + SqrtEpsi * SqrtEps(j) * x * (x - 2.)
		enddo
	enddo
end function

real*8 function LJScoreBond(Pos, Radius, SqrtEps, Cutoff, Bond23, Bond4, &
	& Scale4, NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma, SqrtEpsi
	real*8, dimension(Dim) :: Posi
	real*8 :: Scale
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	LJScoreBond = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				Scale = Scale4
			else
				Scale = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			LJScoreBond = LJScoreBond + Scale * SqrtEpsi * SqrtEps(j) * x * (x - 2.)
		enddo
	enddo
end function

real*8 function LJScoreResBond(Pos, Radius, SqrtEps, Cutoff, ResNum, AtomResNum, &
	& Bond23, Bond4, Scale4, NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim, NBond23, NBond4
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4
	integer :: i, j
	real*8 :: x, CutSq, Radiusi, Sigma, SqrtEpsi, Scale
	real*8, dimension(Dim) :: Posi
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	logical :: IsResi, IsResj
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	LJScoreResBond = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		IsResi = (AtomResNum(i)==ResNum)
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				Scale = Scale4
			else
				Scale = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			IsResj = (AtomResNum(j)==ResNum)
			if (.not. (IsResi .or. IsResj)) cycle
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			Sigma = Radiusi + Radius(j)
			x = Sigma * Sigma / x
			x = x * x * x
			LJScoreResBond = LJScoreResBond + Scale * SqrtEpsi * SqrtEps(j) * x * (x - 2.)
		enddo
	enddo
end function


!======== ENERGIES (LJ + CHARGE) ========

real*8 function Energy(Pos, Radius, SqrtEps, Charge, LJCutoff, ChargeCutoff, CoulFact, NAtom, Dim)
	implicit none
	real*8, intent(in) :: LJCutoff, ChargeCutoff, CoulFact
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps, Charge
	integer :: i, j
	real*8 :: x, DistSq, LJCutSq, ChargeCutSq, Radiusi, Sigma, SqrtEpsi, Chargei
	real*8 :: ELJ, EC
	real*8, dimension(Dim) :: Posi
	ELJ = 0.
	EC = 0.
	LJCutSq = LJCutoff * LJCutoff
	ChargeCutSq = ChargeCutoff * ChargeCutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		Chargei = Charge(i)
		do j = i+1, NAtom
			DistSq = sum((Posi - Pos(j,:))**2)
			if (DistSq < LJCutSq) then
				Sigma = Radiusi + Radius(j)
				x = Sigma * Sigma / DistSq
				x = x * x * x
				ELJ = ELJ + SqrtEpsi * SqrtEps(j) * x * (x - 2.)
			endif
			if (DistSq < ChargeCutSq) then
				EC = EC + Chargei * Charge(j) / sqrt(DistSq)
			endif
		enddo
	enddo
	Energy = ELJ + CoulFact * EC
end function

real*8 function EnergyRes(Pos, Radius, SqrtEps, Charge, LJCutoff, ChargeCutoff, CoulFact, &
	& ResNum, AtomResNum, NAtom, Dim)
	implicit none
	real*8, intent(in) :: LJCutoff, ChargeCutoff, CoulFact
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps, Charge
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer :: i, j
	real*8 :: x, DistSq, LJCutSq, ChargeCutSq, Radiusi, Sigma, SqrtEpsi, Chargei
	real*8 :: ELJ, EC
	real*8, dimension(Dim) :: Posi
	ELJ = 0.
	EC = 0.
	LJCutSq = LJCutoff * LJCutoff
	ChargeCutSq = ChargeCutoff * ChargeCutoff
	do i = 1, NAtom
		if (AtomResNum(i) /= ResNum) cycle
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		Chargei = Charge(i)
		do j = 1, NAtom
			if (AtomResNum(j)==ResNum .and. j <= i) cycle
			DistSq = sum((Posi - Pos(j,:))**2)
			if (DistSq < LJCutSq) then
				Sigma = Radiusi + Radius(j)
				x = Sigma * Sigma / DistSq
				x = x * x * x
				ELJ = ELJ + SqrtEpsi * SqrtEps(j) * x * (x - 2.)
			endif
			if (DistSq < ChargeCutSq) then
				EC = EC + Chargei * Charge(j) / sqrt(DistSq)
			endif
		enddo
	enddo
	EnergyRes = ELJ + CoulFact * EC
end function

real*8 function EnergyBond(Pos, Radius, SqrtEps, Charge, LJCutoff, ChargeCutoff, CoulFact, &
	& Bond23, Bond4, Scale4LJ, Scale4Charge, NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: LJCutoff, ChargeCutoff, CoulFact
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps, Charge
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4LJ, Scale4Charge
	integer :: i, j
	real*8 :: x, DistSq, LJCutSq, ChargeCutSq, Radiusi, Sigma, SqrtEpsi, Chargei
	real*8 :: ELJ, EC
	real*8, dimension(Dim) :: Posi
	real*8 :: ScaleLJ, ScaleCharge
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	ELJ = 0.
	EC = 0.
	LJCutSq = LJCutoff * LJCutoff
	ChargeCutSq = ChargeCutoff * ChargeCutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		Chargei = Charge(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				ScaleLJ = Scale4LJ
				ScaleCharge = Scale4Charge
			else
				ScaleLJ = 1.
				ScaleCharge = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			DistSq = sum((Posi - Pos(j,:))**2)
			if (DistSq < LJCutSq) then
				Sigma = Radiusi + Radius(j)
				x = Sigma * Sigma / DistSq
				x = x * x * x
				ELJ = ELJ + ScaleLJ * SqrtEpsi * SqrtEps(j) * x * (x - 2.)
			endif
			if (DistSq < ChargeCutSq) then
				EC = EC + ScaleCharge * Chargei * Charge(j) / sqrt(DistSq)
			endif
		enddo
	enddo
	EnergyBond = ELJ + CoulFact * EC
end function

real*8 function EnergyResBond(Pos, Radius, SqrtEps, Charge, LJCutoff, ChargeCutoff, CoulFact, &
	& ResNum, AtomResNum, Bond23, Bond4, Scale4LJ, Scale4Charge, NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: LJCutoff, ChargeCutoff, CoulFact
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Radius, SqrtEps, Charge
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4LJ, Scale4Charge
	integer :: i, j
	real*8 :: x, DistSq, LJCutSq, ChargeCutSq, Radiusi, Sigma, SqrtEpsi, Chargei
	real*8 :: ELJ, EC
	real*8, dimension(Dim) :: Posi
	real*8 :: ScaleLJ, ScaleCharge
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	logical :: IsResi, IsResj
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	ELJ = 0.
	EC = 0.
	LJCutSq = LJCutoff * LJCutoff
	ChargeCutSq = ChargeCutoff * ChargeCutoff
	do i = 1, NAtom
		IsResi = (AtomResNum(i) == ResNum)
		Posi = Pos(i,:)
		Radiusi = Radius(i)
		SqrtEpsi = SqrtEps(i)
		Chargei = Charge(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				ScaleLJ = Scale4LJ
				ScaleCharge = Scale4Charge
			else
				ScaleLJ = 1.
				ScaleCharge = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			IsResj = (AtomResNum(j)==ResNum)
			if (.not. (IsResi .or. IsResj)) cycle
			DistSq = sum((Posi - Pos(j,:))**2)
			if (DistSq < LJCutSq) then
				Sigma = Radiusi + Radius(j)
				x = Sigma * Sigma / DistSq
				x = x * x * x
				ELJ = ELJ + ScaleLJ * SqrtEpsi * SqrtEps(j) * x * (x - 2.)
			endif
			if (DistSq < ChargeCutSq) then
				EC = EC + ScaleCharge * Chargei * Charge(j) / sqrt(DistSq)
			endif
		enddo
	enddo
	EnergyResBond = ELJ + CoulFact * EC
end function



!======== CHARGE ENERGIES ========

real*8 function ChargeScore(Pos, Charge, Cutoff, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Charge
	integer :: i, j
	real*8 :: x, CutSq, Chargei
	real*8, dimension(Dim) :: Posi
	ChargeScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Chargei = Charge(i)
		do j = i+1, NAtom
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			ChargeScore = ChargeScore + Chargei * Charge(j) / sqrt(x)
		enddo
	enddo
end function

real*8 function ChargeScoreRes(Pos, Charge, Cutoff, ResNum, AtomResNum, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Charge
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer :: i, j
	real*8 :: x, CutSq, Chargei
	real*8, dimension(Dim) :: Posi
	ChargeScoreRes = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		if (AtomResNum(i) /= ResNum) cycle
		Posi = Pos(i,:)
		Chargei = Charge(i)
		do j = 1, NAtom
			if (AtomResNum(j)==ResNum .and. j <= i) cycle
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			ChargeScoreRes = ChargeScoreRes + Chargei * Charge(j) / sqrt(x)
		enddo
	enddo
end function

real*8 function ChargeScoreBond(Pos, Charge, Cutoff, Bond23, Bond4, Scale4, &
	& NBond23, NBond4, NAtom, Dim)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Charge
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4
	integer :: i, j
	real*8 :: x, CutSq, Chargei
	real*8, dimension(Dim) :: Posi
	real*8 :: Scale
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	ChargeScoreBond = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		Posi = Pos(i,:)
		Chargei = Charge(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				Scale = Scale4
			else
				Scale = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			ChargeScoreBond = ChargeScoreBond + Scale * Chargei * Charge(j) / sqrt(x)
		enddo
	enddo
end function

real*8 function ChargeScoreResBond(Pos, Charge, Cutoff, ResNum, AtomResNum, &
	& Bond23, Bond4, Scale4, NAtom, Dim, NBond23, NBond4)
	implicit none
	real*8, intent(in) :: Cutoff
	integer, intent(in) :: ResNum
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Charge
	integer, dimension(NAtom), intent(in) :: AtomResNum
	integer, intent(in) :: NBond23, NBond4
	integer, dimension(NBond23), intent(in) :: Bond23
	integer, dimension(NBond4), intent(in) :: Bond4
	real*8, intent(in) :: Scale4
	integer :: i, j
	real*8 :: x, CutSq, Chargei
	real*8, dimension(Dim) :: Posi
	real*8 :: Scale
	integer :: CurID, NextID23, Ind23, NextID4, Ind4
	logical :: IsResi, IsResj
	external InitBondIter
	CurID = -1
	call InitBondIter(NextID23, Ind23, Bond23, NAtom, NBond23)
	call InitBondIter(NextID4, Ind4, Bond4, NAtom, NBond4)
	ChargeScoreResBond = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NAtom
		IsResi = (AtomResNum(i) == ResNum)
		Posi = Pos(i,:)
		Chargei = Charge(i)
		do j = i+1, NAtom
			CurID = CurID + 1
			if (CurID == NextID4) then
				if (Ind4 < NBond4) then
					Ind4 = Ind4 + 1
					NextID4 = Bond4(Ind4)
				endif
				Scale = Scale4
			else
				Scale = 1.
			endif
			if (CurID == NextID23) then
				if (Ind23 < NBond23) then
					Ind23 = Ind23 + 1
					NextID23 = Bond23(Ind23)
				endif
				cycle
			endif
			IsResj = (AtomResNum(j)==ResNum)
			if (.not. (IsResi .or. IsResj)) cycle
			x = sum((Posi - Pos(j,:))**2)
			if (x > CutSq) cycle
			ChargeScoreResBond = ChargeScoreResBond + Scale * Chargei * Charge(j) / sqrt(x)
		enddo
	enddo
end function



!======== HYDROGEN BOND SCORES ========

real*8 function HBondScore(Pos, HInd, AccInd, A, B, Cutoff, MinCO, NAtom, Dim, NH, NAcc)
	implicit none
	integer, intent(in) :: NAtom, Dim, NH, NAcc, MinCO
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NH), intent(in) :: HInd
	integer, dimension(NAcc), intent(in) :: AccInd 
	real*8, intent(in) :: A, B, Cutoff
	integer :: i, j
	real*8 :: x, y, E, CutSq
	real*8, dimension(Dim) :: HPos
	HBondScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NH
	    if (HInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		do j = 1, NAcc
			if (abs(i-j) < MinCO) cycle
		    if (AccInd(j) < 0) cycle
			x = sum((HPos - Pos(AccInd(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / x
			y = x * x * x
			E = (A * y - B * x * x) * y
			HBondScore = HBondScore + E
		enddo
	enddo
end function

real*8 function HBondScoreInd(Pos, HInd, AccInd, A, B, Cutoff, MinCO, ResInd, NAtom, Dim, NH, NAcc, NResInd)
	implicit none
	integer, intent(in) :: NAtom, Dim, NH, NAcc, MinCO, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NH), intent(in) :: HInd
	integer, dimension(NAcc), intent(in) :: AccInd 
	real*8, intent(in) :: A, B, Cutoff
	integer, dimension(NResInd), intent(in) :: ResInd
	integer :: m, n, i, j
	real*8 :: x, y, E, CutSq
	real*8, dimension(Dim) :: HPos
	HBondScoreInd = 0.
	CutSq = Cutoff * Cutoff
	do m = 1, NResInd
		i = ResInd(m) + 1
	    if (HInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		do n = 1, NResInd
			j = ResInd(n) + 1
			if (abs(i-j) < MinCO) cycle
		    if (AccInd(j) < 0) cycle
			x = sum((HPos - Pos(AccInd(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / x
			y = x * x * x
			E = (A * y - B * x * x) * y
			HBondScoreInd = HBondScoreInd + E
		enddo
	enddo
end function

real*8 function HBondScore2(Pos1, Pos2, HInd1, AccInd1, HInd2, AccInd2, &
	& A, B, Cutoff, NAtom1, NAtom2, Dim, NH1, NAcc1, NH2, NAcc2)
	implicit none
	integer, intent(in) :: NAtom1, NAtom2, Dim, NH1, NAcc1, NH2, NAcc2
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	integer, dimension(NH1), intent(in) :: HInd1
	integer, dimension(NAcc1), intent(in) :: AccInd1
	integer, dimension(NH2), intent(in) :: HInd2
	integer, dimension(NAcc2), intent(in) :: AccInd2
	real*8, intent(in) :: A, B, Cutoff
	integer :: i, j
	real*8 :: x, y, E, CutSq
	real*8, dimension(Dim) :: HPos
	HBondScore2 = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NH1
	    if (HInd1(i) < 0) cycle
		HPos = Pos1(HInd1(i)+1,:)
		do j = 1, NAcc2
		    if (AccInd2(j) < 0) cycle
			x = sum((HPos - Pos2(AccInd2(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / x
			y = x * x * x
			E = (A * y - B * x * x) * y
			HBondScore2 = HBondScore2 + E
		enddo
	enddo
	do i = 1, NH2
		if (HInd2(i) < 0) cycle
		HPos = Pos2(HInd2(i)+1,:)
		do j = 1, NAcc1
			if (AccInd1(j) < 0) cycle
			x = sum((HPos - Pos1(AccInd1(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / x
			y = x * x * x
			E = (A * y - B * x * x) * y
			HBondScore2 = HBondScore2 + E
		enddo
	enddo
end function


real*8 function HBondLJScore(Pos, NInd, HInd, OInd, A, B, Cutoff, CosPower, MinCO, NAtom, Dim, NRes)
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, MinCO
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NRes), intent(in) :: NInd, HInd, OInd
	real*8, intent(in) :: A, B, Cutoff
	integer, intent(in) :: CosPower
	integer :: i, j
	real*8 :: x, y, z, E, CutSq
	real*8, dimension(Dim) :: HPos, NPos, OPos
	HBondLJScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NRes
	    if (HInd(i) < 0 .or. NInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		NPos = Pos(NInd(i)+1,:)
		do j = 1, NRes
			if (abs(i-j) < MinCO) cycle
		    if (OInd(j) < 0) cycle
			OPos = Pos(OInd(j)+1,:)
			x = sum((HPos - OPos)**2)
			if (x > CutSq) cycle
			x = 1. / x
			y = x * x * x
			if (CosPower > 0) then
				z = dot_product(HPos-NPos, OPos-HPos)
				if (z < 0.) then
					E = 0.
				else
					z = z*z * x / sum((HPos - NPos)**2)
					if (CosPower /= 2) z = sqrt(z)**CosPower
					E = (A * y - B * x * x) * y * z
				endif
			else
				E = (A * y - B * x * x) * y
			endif
			HBondLJScore = HBondLJScore + E
		enddo
	enddo
end function

real*8 function HBondLJScoreInd(Pos, NInd, HInd, OInd, A, B, Cutoff, CosPower, MinCO, ResInd, NAtom, Dim, NRes, NResInd)
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, MinCO, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NRes), intent(in) :: NInd, HInd, OInd
	real*8, intent(in) :: A, B, Cutoff
	integer, intent(in) :: CosPower
	integer, dimension(NResInd), intent(in) :: ResInd
	integer :: i, j, m, n
	real*8 :: x, y, z, E, CutSq
	real*8, dimension(Dim) :: HPos, NPos, OPos
	HBondLJScoreInd = 0.
	CutSq = Cutoff * Cutoff
	do m = 1, NResInd
		i = ResInd(m) + 1
	    if (HInd(i) < 0 .or. NInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		NPos = Pos(NInd(i)+1,:)
		do n = 1, NResInd
			j = ResInd(n) + 1
			if (abs(i-j) < MinCO) cycle
		    if (OInd(j) < 0) cycle
			OPos = Pos(OInd(j)+1,:)
			x = sum((HPos - OPos)**2)
			if (x > CutSq) cycle
			x = 1. / x
			y = x * x * x
			if (CosPower > 0) then
				z = dot_product(HPos-NPos, OPos-HPos)
				if (z < 0.) then
					E = 0.
				else
					z = z*z * x / sum((HPos - NPos)**2)
					if (CosPower /= 2) z = sqrt(z)**CosPower
					E = (A * y - B * x * x) * y * z
				endif
			else
				E = (A * y - B * x * x) * y
			endif
			HBondLJScoreInd = HBondLJScoreInd + E
		enddo
	enddo
end function


real*8 function HBondDipoleScore(Pos, NInd, HInd, CInd, OInd, Cutoff, MinCO, NAtom, Dim, NN, NH, NC, NO)
	implicit none
	integer, intent(in) :: NAtom, Dim, NN, NH, NC, NO, MinCO
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NN), intent(in) :: NInd	
	integer, dimension(NH), intent(in) :: HInd
	integer, dimension(NC), intent(in) :: CInd
	integer, dimension(NO), intent(in) :: OInd 
	real*8, intent(in) :: Cutoff
	integer :: i, j
	real*8 :: x, E, CutSq
	real*8, dimension(Dim) :: HPos, NHVec
	HBondDipoleScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NH
		if (HInd(i) < 0 .or. NInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		NHVec = HPos - Pos(NInd(i)+1,:)
		do j = 1, NO
			if (abs(i-j) < MinCO) cycle
			if (OInd(i) < 0 .or. CInd(i) < 0) cycle
			x = sum((HPos - Pos(OInd(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / sqrt(x)
			x = x * x * x
			E = dot_product(NHVec, Pos(OInd(j)+1,:) - Pos(CInd(j)+1,:)) * x
			HBondDipoleScore = HBondDipoleScore + E
		enddo
	enddo
end function

real*8 function HBondDipoleScoreInd(Pos, NInd, HInd, CInd, OInd, Cutoff, MinCO, ResInd, NAtom, Dim, NN, NH, NC, NO, NResInd)
	implicit none
	integer, intent(in) :: NAtom, Dim, NN, NH, NC, NO, MinCO, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NN), intent(in) :: NInd	
	integer, dimension(NH), intent(in) :: HInd
	integer, dimension(NC), intent(in) :: CInd
	integer, dimension(NO), intent(in) :: OInd 
	integer, dimension(NResInd), intent(in) :: ResInd
	real*8, intent(in) :: Cutoff
	integer :: i, j, m, n
	real*8 :: x, E, CutSq
	real*8, dimension(Dim) :: HPos, NHVec
	HBondDipoleScoreInd = 0.
	CutSq = Cutoff * Cutoff
	do m = 1, NResInd
		i = ResInd(m) + 1
		if (HInd(i) < 0 .or. NInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		NHVec = HPos - Pos(NInd(i)+1,:)
		do n = 1, NResInd
			j = ResInd(n) + 1
			if (abs(i-j) < MinCO) cycle
			if (OInd(i) < 0 .or. CInd(i) < 0) cycle
			x = sum((HPos - Pos(OInd(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / sqrt(x)
			x = x * x * x
			E = dot_product(NHVec, Pos(OInd(j)+1,:) - Pos(CInd(j)+1,:)) * x
			HBondDipoleScoreInd = HBondDipoleScoreInd + E
		enddo
	enddo
end function

real*8 function HBondDipoleScore2(Pos1, Pos2, NInd1, HInd1, CInd1, OInd1, &
	& NInd2, HInd2, CInd2, OInd2, Cutoff, NAtom1, NAtom2, Dim, &
	& NN1, NH1, NC1, NO1, NN2, NH2, NC2, NO2)
	implicit none
	integer, intent(in) :: NAtom1, NAtom2, Dim, NN1, NH1, NC1, NO1, NN2, NH2, NC2, NO2
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	integer, dimension(NN1), intent(in) :: NInd1
	integer, dimension(NH1), intent(in) :: HInd1
	integer, dimension(NC1), intent(in) :: CInd1
	integer, dimension(NO1), intent(in) :: OInd1
	integer, dimension(NN2), intent(in) :: NInd2
	integer, dimension(NH2), intent(in) :: HInd2
	integer, dimension(NC2), intent(in) :: CInd2
	integer, dimension(NO2), intent(in) :: OInd2	
	real*8, intent(in) :: Cutoff
	integer :: i, j
	real*8 :: x, E, CutSq
	real*8, dimension(Dim) :: HPos, NHVec
	HBondDipoleScore2 = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NH1
		if (HInd1(i) < 0 .or. NInd1(i) < 0) cycle
		HPos = Pos1(HInd1(i)+1,:)
		NHVec = HPos - Pos1(NInd1(i)+1,:)
		do j = 1, NO2
			if (OInd2(i) < 0 .or. CInd2(i) < 0) cycle
			x = sum((HPos - Pos2(OInd2(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / sqrt(x)
			x = x * x * x
			E = dot_product(NHVec, Pos2(OInd2(j)+1,:) - Pos2(CInd2(j)+1,:)) * x
			HBondDipoleScore2 = HBondDipoleScore2 + E
		enddo
	enddo
	do i = 1, NH2
		if (HInd2(i) < 0 .or. NInd2(i) < 0) cycle
		HPos = Pos2(HInd2(i)+1,:)
		NHVec = HPos - Pos2(NInd2(i)+1,:)
		do j = 1, NO1
			if (OInd1(i) < 0 .or. CInd1(i) < 0) cycle
			x = sum((HPos - Pos1(OInd1(j)+1,:))**2)
			if (x > CutSq) cycle
			x = 1. / sqrt(x)
			x = x * x * x
			E = dot_product(NHVec, Pos1(OInd1(j)+1,:) - Pos1(CInd1(j)+1,:)) * x
			HBondDipoleScore2 = HBondDipoleScore2 + E
		enddo
	enddo
end function



real*8 function HBondChargeScore(Pos, NInd, HInd, CInd, OInd, Cutoff, &
	& CosPower1, CosPower2, Coef, MinCO, NAtom, Dim, NRes)
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, MinCO
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NRes), intent(in) :: NInd, HInd, CInd, OInd	
	real*8, intent(in) :: Cutoff, Coef
	integer, intent(in) :: CosPower1, CosPower2
	integer :: i, j
	real*8 :: dOH, dON, dCH, dCN, dHN, dCO, E, CutSq, Cos1, Cos2
	real*8, dimension(Dim) :: HPos, NPos, CPos, OPos
	HBondChargeScore = 0.
	CutSq = Cutoff * Cutoff
	do i = 1, NRes
		if (HInd(i) < 0 .or. NInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		NPos = Pos(NInd(i)+1,:)
		do j = 1, NRes
			if (abs(i-j) < MinCO) cycle
			if (OInd(j) < 0 .or. CInd(j) < 0) cycle
			OPos = Pos(OInd(j)+1,:)
			CPos = Pos(CInd(j)+1,:)
			dOH = sum((HPos - OPos)**2)
			if (dOH > CutSq) cycle
			dOH = 1. / sqrt(dOH)
			dON = 1. / sqrt(sum((NPos - OPos)**2))
			dCH = 1. / sqrt(sum((HPos - CPos)**2))
			dCN = 1. / sqrt(sum((NPos - CPos)**2))
			dHN = 1. / sqrt(sum((NPos - HPos)**2))
			dCO = 1. / sqrt(sum((OPos - CPos)**2))
			E = Coef*(dON + dCH - dCN) - dOH
			Cos1 = max(0., dot_product(HPos - OPos, NPos - HPos) * dOH * dHN)
			Cos2 = max(0., dot_product(HPos - OPos, OPos - CPos) * dOH * dCO)
			HBondChargeScore = HBondChargeScore + E * (Cos1**CosPower1) * (Cos2**CosPower2)
		enddo
	enddo
end function

real*8 function HBondChargeScoreInd(Pos, NInd, HInd, CInd, OInd, Cutoff, &
	& CosPower1, CosPower2, Coef, MinCO, ResInd, NAtom, Dim, NRes, NResInd)
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, MinCO, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NRes), intent(in) :: NInd, HInd, CInd, OInd
	integer, dimension(NResInd), intent(in) :: ResInd
	real*8, intent(in) :: Cutoff, Coef
	integer, intent(in) :: CosPower1, CosPower2
	integer :: i, j, m, n
	real*8 :: dOH, dON, dCH, dCN, dHN, dCO, E, CutSq, Cos1, Cos2
	real*8, dimension(Dim) :: HPos, NPos, CPos, OPos
	HBondChargeScoreInd = 0.
	CutSq = Cutoff * Cutoff
	do m = 1, NResInd
		i = ResInd(m) + 1
		if (HInd(i) < 0 .or. NInd(i) < 0) cycle
		HPos = Pos(HInd(i)+1,:)
		NPos = Pos(NInd(i)+1,:)
		do n = 1, NResInd
			j = ResInd(n) + 1
			if (abs(i-j) < MinCO) cycle
			if (OInd(j) < 0 .or. CInd(j) < 0) cycle
			OPos = Pos(OInd(j)+1,:)
			CPos = Pos(CInd(j)+1,:)
			dOH = sum((HPos - OPos)**2)
			if (dOH > CutSq) cycle
			dOH = 1. / sqrt(dOH)
			dON = 1. / sqrt(sum((NPos - OPos)**2))
			dCH = 1. / sqrt(sum((HPos - CPos)**2))
			dCN = 1. / sqrt(sum((NPos - CPos)**2))
			dHN = 1. / sqrt(sum((NPos - HPos)**2))
			dCO = 1. / sqrt(sum((OPos - CPos)**2))
			E = Coef*(dON + dCH - dCN) - dOH
			Cos1 = max(0., dot_product(HPos - OPos, NPos - HPos) * dOH * dHN)
			Cos2 = max(0., dot_product(HPos - OPos, OPos - CPos) * dOH * dCO)
			HBondChargeScoreInd = HBondChargeScoreInd + E * (Cos1**CosPower1) * (Cos2**CosPower2)
		enddo
	enddo
end function


!======== RG CALCULATIONS ========

real*8 function Rg(Pos, NAtom, Dim)
	implicit none
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(Dim) :: Center
	integer :: i
	Center = sum(Pos, 1) / real(NAtom)
	Rg = 0.
	do i = 1, NAtom
		Rg = Rg + sum((Pos(i,:) - Center)**2)
	enddo
	Rg = Rg / float(NAtom)
	Rg = sqrt(Rg)
end function

real*8 function RgWeights(Pos, Weights, NAtom, Dim)
	implicit none
	integer, intent(in) :: NAtom, Dim
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	real*8, dimension(NAtom), intent(in) :: Weights
	real*8, dimension(Dim) :: Center
	integer :: i
	Center = sum(Pos, 1) / real(NAtom)
	RgWeights = 0.
	do i = 1, NAtom
		RgWeights = RgWeights + Weights(i) * sum((Pos(i,:) - Center)**2)
	enddo
	RgWeights = RgWeights / sum(Weights)
	RgWeights = sqrt(RgWeights)
end function

real*8 function RgAtomInd(Pos, AtomInd, NAtom, Dim, NInd)
	implicit none
	integer, intent(in) :: NAtom, Dim, NInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NInd), intent(in) :: AtomInd
	real*8, dimension(Dim) :: Center
	integer :: i
	Center = 0.
	RgAtomInd = 0.
	do i = 1, NInd
		Center = Center + Pos(AtomInd(i) + 1, :)
	enddo
	Center = Center / real(NInd)
	do i = 1, NInd
		RgAtomInd = RgAtomInd + sum((Pos(AtomInd(i) + 1,:) - Center)**2)
	enddo
	RgAtomInd = RgAtomInd / real(NInd)
	RgAtomInd = sqrt(RgAtomInd)
end function


real*8 function Rg2(Pos1, Pos2, NAtom1, NAtom2, Dim)
	implicit none
	integer, intent(in) :: NAtom1, NAtom2, Dim
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	real*8, dimension(Dim) :: Center
	integer :: i
	Center = (sum(Pos1, 1) + sum(Pos2, 1)) / real(NAtom1 + NAtom2)
	Rg2 = 0.
	do i = 1, NAtom1
		Rg2 = Rg2 + sum((Pos1(i,:) - Center)**2)
	enddo
	do i = 1, NAtom2
		Rg2 = Rg2 + sum((Pos2(i,:) - Center)**2)
	enddo
	Rg2 = Rg2 / real(NAtom1 + NAtom2)
	Rg2 = sqrt(Rg2)
end function

real*8 function RgWeights2(Pos1, Pos2, Weights1, Weights2, NAtom1, NAtom2, Dim)
	implicit none
	integer, intent(in) :: NAtom1, NAtom2, Dim
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	real*8, dimension(NAtom1), intent(in) :: Weights1
	real*8, dimension(NAtom2), intent(in) :: Weights2
	real*8, dimension(Dim) :: Center
	integer :: i
	Center = (sum(Pos1, 1) + sum(Pos2, 1)) / real(NAtom1 + NAtom2)
	RgWeights2 = 0.
	do i = 1, NAtom1
		RgWeights2 = RgWeights2 + Weights1(i) * sum((Pos1(i,:) - Center)**2)
	enddo
	do i = 1, NAtom2
		RgWeights2 = RgWeights2 + Weights2(i) * sum((Pos2(i,:) - Center)**2)
	enddo
	RgWeights2 = RgWeights2 / (sum(Weights1) + sum(Weights2))
	RgWeights2 = sqrt(RgWeights2)
end function

real*8 function RgAtomInd2(Pos1, Pos2, AtomInd1, AtomInd2, NAtom1, NAtom2, Dim, NInd1, NInd2)
	implicit none
	integer, intent(in) :: NAtom1, NAtom2, Dim, NInd1, NInd2
	real*8, dimension(NAtom1,Dim), intent(in) :: Pos1
	real*8, dimension(NAtom2,Dim), intent(in) :: Pos2
	integer, dimension(NInd1), intent(in) :: AtomInd1
	integer, dimension(NInd2), intent(in) :: AtomInd2
	real*8, dimension(Dim) :: Center
	integer :: i
	Center = 0.
	RgAtomInd2 = 0.
	do i = 1, NInd1
		Center = Center + Pos1(AtomInd1(i) + 1, :)
	enddo
	do i = 1, NInd2
		Center = Center + Pos2(AtomInd2(i) + 1, :)
	enddo
	Center = Center / real(NInd1 + NInd2)
	do i = 1, NInd1
		RgAtomInd2 = RgAtomInd2 + sum((Pos1(AtomInd1(i) + 1,:) - Center)**2)
	enddo
	do i = 1, NInd2
		RgAtomInd2 = RgAtomInd2 + sum((Pos2(AtomInd2(i) + 1,:) - Center)**2)
	enddo
	RgAtomInd2 = RgAtomInd2 / real(NInd1 + NInd2)
	RgAtomInd2 = sqrt(RgAtomInd2)
end function


!======== BOND TESTERS AND MAPS ========

integer function IsBond(Pos1, Pos2, Pos3, DistCut, AngCut, Dim)
	!here we are testing for a bond from 1 to 2, where 3 is already bonded to 2
	implicit none
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: Pos1, Pos2, Pos3
	real*8, intent(in) :: DistCut, AngCut
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: RadPerDeg = pi / 180.D0
	real*8 :: DistSq, CosAng
	real*8, dimension(Dim) :: Dist32, Dist21
	IsBond = 0
	!first check distance
	Dist21 = Pos1 - Pos2 
	DistSq = dot_product(Dist21,Dist21)
	if (DistSq > DistCut*DistCut) return
	!now check angle
	Dist32 = Pos2 - Pos3
	DistSq = DistSq * dot_product(Dist32,Dist32)
	DistSq = max(DistSq, tiny(DistSq))
	CosAng = dot_product(Dist32, Dist21) / sqrt(DistSq)
	if (CosAng > cos(AngCut * RadPerDeg)) return
	IsBond = 1
end function


integer function IsHBond(NPos, HPos, CPos, OPos, ECutoff, Dim)
	!DSSP-type hbond assignment
	implicit none
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: NPos, HPos, CPos, OPos
	real*8, intent(in) :: ECutoff
	real*8 :: dOH, dON, dCH, dCN
	dOH = sqrt(sum((HPos - OPos)**2))
	dON = sqrt(sum((NPos - OPos)**2))
	dCH = sqrt(sum((HPos - CPos)**2))
	dCN = sqrt(sum((NPos - CPos)**2))
	if (dOH > dON .or. dOH > dCH .or. dOH > dCN) then
		IsHBond = 0
	elseif (1./dON + 1./dCH - 1./dCN - 1./dOH < ECutoff) then
		IsHBond = 1
	else
		IsHBond = 0
	endif
end function



subroutine BondContactMap(Pos, Ind1, Ind2, Ind3, ResInd, DistCut, AngCut, CM, NAtom, Dim, NRes, NResInd)
	!here we are testing for a bond from 1 to 2, where 3 is already bonded to 2
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NRes), intent(in) :: Ind1, Ind2, Ind3
	integer, dimension(NResInd), intent(in) :: ResInd
	real*8, intent(in) :: DistCut, AngCut
	integer, dimension(NRes,NRes), intent(out) :: CM
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: RadPerDeg = pi / 180.D0
	integer :: i, j, r1, r2, a1, a2, a3
	real*8 :: DistCutSq, DistSq, CosAngCut, CosAng
	real*8, dimension(Dim) :: Dist32, Dist21
	CM = 0.
	DistCutSq = DistCut*DistCut
	CosAngCut = cos(AngCut * RadPerDeg)
	do i = 1, NResInd
		r2 = ResInd(i) + 1
		a2 = Ind2(r2) + 1
		a3 = Ind3(r2) + 1
		if (a2 <= 0 .or. a3 <= 0) cycle
		Dist32 = Pos(a2,:) - Pos(a3,:)
		do j = 1, NResInd
			r1 = ResInd(j) + 1
			if (r1==r2) cycle
			a1 = Ind1(r1) + 1
			if (a1 <= 0) cycle
			!first check distance
			Dist21 = Pos(a1,:) - Pos(a2,:) 
			DistSq = dot_product(Dist21,Dist21)
			if (DistSq > DistCutSq) cycle
			!now check angle
			DistSq = DistSq * dot_product(Dist32,Dist32)
			DistSq = max(DistSq, tiny(DistSq))
			CosAng = dot_product(Dist32, Dist21) / sqrt(DistSq)
			if (CosAng < CosAngCut) cycle
			CM(r1,r2) = 1
		enddo
	enddo
end subroutine


subroutine HBondContactMap(Pos, NInd, HInd, CInd, OInd, ResInd, ECutoff, CM, NAtom, Dim, NRes, NResInd)
	!here we are testing for a bond from 1 to 2, where 3 is already bonded to 2
	!this uses a dssp approach
	implicit none
	integer, intent(in) :: NAtom, Dim, NRes, NResInd
	real*8, dimension(NAtom,Dim), intent(in) :: Pos
	integer, dimension(NRes), intent(in) :: NInd, HInd, CInd, OInd
	integer, dimension(NResInd), intent(in) :: ResInd
	real*8, intent(in) :: ECutoff
	integer, dimension(NRes,NRes), intent(out) :: CM
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: RadPerDeg = pi / 180.D0
	integer :: i, j, rj, ri
	real*8, dimension(Dim) :: NPos, HPos, CPos, OPos
	integer, external :: IsHBond
	CM = 0.
	do i = 1, NResInd
		ri = ResInd(i) + 1
		if (NInd(ri) < 0 .or. HInd(ri) < 0) cycle
		NPos = Pos(NInd(ri)+1,:)
		HPos = Pos(HInd(ri)+1,:)
		do j = 1, NResInd
			rj = ResInd(j) + 1
			if (ri==rj) cycle
			if (CInd(rj) < 0 .or. OInd(rj) < 0) cycle
			CPos = Pos(CInd(rj)+1,:)
			OPos = Pos(OInd(rj)+1,:)
			CM(ri,rj) = real(IsHBond(NPos, HPos, CPos, OPos, ECutoff, Dim))
		enddo
	enddo
end subroutine
	 

!======== MISC FUNCTIONS ========

integer function GetRandomInd(Prob, RanNum, N)
	integer, intent(in) :: N
	real*8, dimension(N), intent(in) :: Prob
	real*8, intent(in) :: RanNum
	real*8 :: CumProb
	CumProb = Prob(1)
	GetRandomInd = 1
	do while (GetRandomInd < N .and. RanNum >= CumProb)
		GetRandomInd = GetRandomInd + 1
		CumProb = CumProb + Prob(GetRandomInd)
	enddo
	!correct for python ordering
	GetRandomInd = GetRandomInd - 1
end function
