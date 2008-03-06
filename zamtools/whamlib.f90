
subroutine free_energy(EKN, BetaK, NBin, NIterBin, NIterAll, FK, NTemp, NFrame)
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: EKN
	real*8, dimension(NTemp), intent(in) :: BetaK
	integer, intent(in) :: NBin, NIterBin, NIterAll
	real*8, dimension(NTemp), intent(out) :: FK
	real*8 :: LogNFrame
	real*8 :: EMin, EMax, EDelta, InvEDelta
	real*8, dimension(:), allocatable :: LogHist, LogP, HistBinVal
	real*8, dimension(:,:), allocatable :: LogPKN
	real*8, dimension(NTemp) :: LogDem, EAvgK
	real*8 :: LogNum, LogDemMax, LogPMax, E
	integer :: i, j, k, k2, n

	!make allocations
	allocate(LogHist(NBin))
	allocate(LogP(NBin))
	allocate(HistBinVal(NBin))
	if (NIterAll > 0) allocate(LogPKN(NTemp, NFrame))

	!create histograms
	EMin = minval(EKN)
	EMax = maxval(EKN)
	EDelta = max((EMax - EMin) / real(NBin), 1.d-100)
	InvEDelta = 1. / EDelta
	LogHist = (100. * tiny(LogHist(1)))
	do k = 1, NTemp
		do n = 1, NFrame
			i = min(int((EKN(k,n) - EMin) * InvEDelta) + 1, NBin)
			LogHist(i) = LogHist(i) + 1.
		enddo
	enddo
	do j = 1, NBin
		HistBinVal(j) = EMin + EDelta * (real(j) - 0.5)
	enddo

	!take logs
	LogHist = log(LogHist)
	LogNFrame = log(real(NFrame))

	!find average energies
	EAvgK = sum(EKN, 2) / real(NFrame)

	!create initial prob and free energy est
	FK = 0.

	!bin iteration to refine FK
	do i = 1, NIterBin
		do k = 1, NTemp
			do j = 1, NBin
				LogNum = LogHist(j) - BetaK(k) * HistBinVal(j)
				LogDem = FK - BetaK * HistBinVal(j)
				LogDemMax = maxval(LogDem)
				LogP(j) = LogNum - log(sum(exp(LogDem-LogDemMax))) - LogDemMax - LogNFrame
			enddo
			LogPMax = maxval(LogP)
			FK(k) = -log(sum(exp(LogP - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
	enddo
	!cleanup
	deallocate(LogHist)
	deallocate(LogP)
	deallocate(HistBinVal)

	!iteration over all data (time consuming; ever so slightly more accurate)
	do i = 1, NIterAll
		do k = 1, NTemp
			do k2 = 1, NTemp
				do n = 1, NFrame
					E = EKN(k2,n)
					LogNum = -BetaK(k) * E
					LogDem = FK - BetaK * E
					LogDemMax = maxval(LogDem)
					LogPKN(k2,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
				enddo
			enddo
			LogPMax = maxval(LogPKN)
			FK(k) = -log(sum(exp(LogPKN - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
	enddo

	!cleanup
	if (allocated(LogPKN)) deallocate(LogPKN)
end subroutine


subroutine log_weight(EKN, BetaK, TargetBeta, FK, LogwKN, NTemp, NFrame)
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: EKN
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, intent(in) :: TargetBeta
	real*8, dimension(NTemp), intent(in) :: FK
	real*8, dimension(NTemp, NFrame), intent(out) :: LogwKN
	real*8 :: LogNFrame
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax
	integer :: k, n
	!calculate individual probabilities
	LogNFrame = log(real(NFrame))
	do k = 1, NTemp
		do n = 1, NFrame
			LogNum = -TargetBeta * EKN(k,n)
			LogDem = FK - BetaK * EKN(k,n)
			LogDemMax = maxval(LogDem)
			LogwKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
		enddo
	enddo
	!shift probabilities
	LogwKN = LogwKN - maxval(LogwKN)	
end subroutine


subroutine log_weight_scaleterm(ETotKN, ECompKN, Scale, BetaK, TargetBeta, FK, LogwKN, NTemp, NFrame)
	!This will calculate the log weights when a component term in the energy function is
	!scaled to a different value than during the run.
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: ETotKN, ECompKN
	real*8, intent(in) :: Scale
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, intent(in) :: TargetBeta
	real*8, dimension(NTemp), intent(in) :: FK
	real*8, dimension(NTemp, NFrame), intent(out) :: LogwKN
	real*8 :: LogNFrame
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax
	integer :: k, n
	!calculate individual probabilities
	LogNFrame = log(real(NFrame))
	do k = 1, NTemp
		do n = 1, NFrame
			LogNum = -TargetBeta * (ETotKN(k,n) + (Scale - 1.) * ECompKN(k,n))
			LogDem = FK - BetaK * ETotKN(k,n)
			LogDemMax = maxval(LogDem)
			LogwKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
		enddo
	enddo
	!shift probabilities
	LogwKN = LogwKN - maxval(LogwKN)	
end subroutine

subroutine log_weight_restraint(ETotKN, ERestKN, BetaK, TargetBeta, FK, LogwKN, NTemp, NFrame)
	!This will calculate the log weights for a restrained simulation when the restraint
	!is removed.  ETotKN is the total energy (incl constraint) and ERestKN is the cons energy.
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: ETotKN, ERestKN
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, intent(in) :: TargetBeta
	real*8, dimension(NTemp), intent(in) :: FK
	real*8, dimension(NTemp, NFrame), intent(out) :: LogwKN
	real*8 :: LogNFrame
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax
	integer :: k, n
	!calculate individual probabilities
	LogNFrame = log(real(NFrame))
	do k = 1, NTemp
		do n = 1, NFrame
			LogNum = -TargetBeta * (ETotKN(k,n) - ERestKN(k,n))
			LogDem = FK - BetaK * ETotKN(k,n)
			LogDemMax = maxval(LogDem)
			LogwKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
		enddo
	enddo
	!shift probabilities
	LogwKN = LogwKN - maxval(LogwKN)	
end subroutine


!======== TWO TERMS IN THE ENERGY FUNCTION ========

subroutine free_energy_2(EKN1, EKN2, WeightK1, WeightK2, BetaK, &
	& NBin, NIterBin, NIterAll, FK, NTemp, NFrame)
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: EKN1, EKN2
	real*8, dimension(NTemp), intent(in) :: WeightK1, WeightK2
	real*8, dimension(NTemp), intent(in) :: BetaK
	integer, intent(in) :: NBin, NIterBin, NIterAll
	real*8, dimension(NTemp), intent(out) :: FK
	real*8 :: LogNFrame
	real*8 :: EMin1, EMax1, EDelta1, InvEDelta1
	real*8 :: EMin2, EMax2, EDelta2, InvEDelta2
	real*8, dimension(:), allocatable :: HistBinVal1, HistBinVal2
	real*8, dimension(:,:), allocatable :: LogHist, LogP
	real*8, dimension(:,:), allocatable :: LogPKN
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax, LogPMax, E1, E2
	integer :: i, j, k, k2, l, n
	!note: we're making explicit duplicate varibles (XX1, XX2, etc) to avoid multiple indices, which are SLOW

	!make allocations
	allocate(LogHist(NBin,NBin))
	allocate(LogP(NBin,NBin))
	allocate(HistBinVal1(NBin), HistBinVal2(NBin))
	if (NIterAll > 0) allocate(LogPKN(NTemp, NFrame))

	!create histograms
	EMin1 = minval(EKN1)
	EMax1 = maxval(EKN1)
	EDelta1 = max((EMax1 - EMin1) / real(NBin), 1.d-100)
	InvEDelta1 = 1. / EDelta1
	EMin2 = minval(EKN2)
	EMax2 = maxval(EKN2)
	EDelta2 = max((EMax2 - EMin2) / real(NBin), 1.d-100)
	InvEDelta2 = 1. / EDelta2

	LogHist = (100. * tiny(LogHist(1,1)))
	do k = 1, NTemp
		do n = 1, NFrame
			i = min(int((EKN1(k,n) - EMin1) * InvEDelta1) + 1, NBin)
			j = min(int((EKN2(k,n) - EMin2) * InvEDelta2) + 1, NBin)
			LogHist(i,j) = LogHist(i,j) + 1.
		enddo
	enddo

	do j = 1, NBin
		HistBinVal1(j) = EMin1 + EDelta1 * (real(j) - 0.5)
		HistBinVal2(j) = EMin2 + EDelta2 * (real(j) - 0.5)
	enddo

	!take logs
	LogHist = log(LogHist)
	LogNFrame = log(real(NFrame))

	!create initial prob and free energy est
	FK = 0.

	!bin iteration to refine FK
	do i = 1, NIterBin
		do k = 1, NTemp
			do j = 1, NBin
				do l = 1, NBin
					LogNum = LogHist(j,l) - BetaK(k) * (WeightK1(k) * HistBinVal1(j) + WeightK2(k) * HistBinVal2(l))
					LogDem = FK - BetaK * (WeightK1 * HistBinVal1(j) + WeightK2 * HistBinVal2(l))
					LogDemMax = maxval(LogDem)
					LogP(j,l) = LogNum - log(sum(exp(LogDem-LogDemMax))) - LogDemMax - LogNFrame
				enddo
			enddo
			LogPMax = maxval(LogP)
			FK(k) = -log(sum(exp(LogP - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
	enddo
	!cleanup
	deallocate(LogHist)
	deallocate(LogP)
	deallocate(HistBinVal1, HistBinVal2)

	!iteration over all data (time consuming; ever so slightly more accurate)
	do i = 1, NIterAll
		do k = 1, NTemp
			do k2 = 1, NTemp
				do n = 1, NFrame
					E1 = EKN1(k2,n)
					E2 = EKN2(k2,n)
					LogNum = -BetaK(k) * (WeightK1(k) * E1 + WeightK2(k) * E2)
					LogDem = FK - BetaK * (WeightK1 * E1 + WeightK2 * E2)
					LogDemMax = maxval(LogDem)
					LogPKN(k2,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
				enddo
			enddo
			LogPMax = maxval(LogPKN)
			FK(k) = -log(sum(exp(LogPKN - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
	enddo

	!cleanup
	if (allocated(LogPKN)) deallocate(LogPKN)
end subroutine


subroutine log_weight_2(EKN1, EKN2, WeightK1, WeightK2, BetaK, TargetWeight1, TargetWeight2, &
	& TargetBeta, FK, LogwKN, NTemp, NFrame)
	!This will calculate the log weights for general weights in the energy function.
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: EKN1, EKN2
	real*8, dimension(NTemp), intent(in) :: WeightK1, WeightK2
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, intent(in) :: TargetWeight1, TargetWeight2, TargetBeta
	real*8, dimension(NTemp), intent(in) :: FK
	real*8, dimension(NTemp, NFrame), intent(out) :: LogwKN
	real*8 :: LogNFrame
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax
	integer :: k, n
	!calculate individual probabilities
	LogNFrame = log(real(NFrame))
	do k = 1, NTemp
		do n = 1, NFrame
			LogNum = -TargetBeta * (EKN1(k,n) * TargetWeight1 + EKN2(k,n) * TargetWeight2)
			LogDem = FK - BetaK * (EKN1(k,n) * WeightK1  + EKN2(k,n) * WeightK2)
			LogDemMax = maxval(LogDem)
			LogwKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
		enddo
	enddo
	!shift probabilities
	LogwKN = LogwKN - maxval(LogwKN)	
end subroutine


!======== THREE TERMS IN THE ENERGY FUNCTION ========

subroutine free_energy_3(EKN1, EKN2, EKN3, WeightK1, WeightK2, WeightK3, &
	& BetaK, NBin, NIterBin, NIterAll, FK, NTemp, NFrame)
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: EKN1, EKN2, EKN3
	real*8, dimension(NTemp), intent(in) :: WeightK1, WeightK2, WeightK3
	real*8, dimension(NTemp), intent(in) :: BetaK
	integer, intent(in) :: NBin, NIterBin, NIterAll
	real*8, dimension(NTemp), intent(out) :: FK
	real*8 :: LogNFrame
	real*8 :: EMin1, EMax1, EDelta1, InvEDelta1
	real*8 :: EMin2, EMax2, EDelta2, InvEDelta2
	real*8 :: EMin3, EMax3, EDelta3, InvEDelta3
	real*8, dimension(:), allocatable :: HistBinVal1, HistBinVal2, HistBinVal3
	real*8, dimension(:,:,:), allocatable :: LogHist, LogP
	real*8, dimension(:,:), allocatable :: LogPKN
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax, LogPMax, E1, E2, E3
	integer :: i1, i2, i3, i, j, k, k2, n
	!note: we're making explicit duplicate varibles (XX1, XX2, etc) to avoid multiple indices, which are SLOW

	!make allocations
	allocate(LogHist(NBin,NBin,NBin))
	allocate(LogP(NBin,NBin,NBin))
	allocate(HistBinVal1(NBin), HistBinVal2(NBin), HistBinVal3(NBin))
	if (NIterAll > 0) allocate(LogPKN(NTemp, NFrame))

	!create histograms
	EMin1 = minval(EKN1)
	EMax1 = maxval(EKN1)
	EDelta1 = max((EMax1 - EMin1) / real(NBin), 1.d-100)
	InvEDelta1 = 1. / EDelta1
	EMin2 = minval(EKN2)
	EMax2 = maxval(EKN2)
	EDelta2 = max((EMax2 - EMin2) / real(NBin), 1.d-100)
	InvEDelta2 = 1. / EDelta2
	EMin3 = minval(EKN2)
	EMax3 = maxval(EKN2)
	EDelta3 = max((EMax3 - EMin3) / real(NBin), 1.d-100)
	InvEDelta3 = 1. / EDelta2

	LogHist = (100. * tiny(LogHist(1,1,1)))
	do k = 1, NTemp
		do n = 1, NFrame
			i1 = min(int((EKN1(k,n) - EMin1) * InvEDelta1) + 1, NBin)
			i2 = min(int((EKN2(k,n) - EMin2) * InvEDelta2) + 1, NBin)
			i3 = min(int((EKN3(k,n) - EMin3) * InvEDelta3) + 1, NBin)
			LogHist(i1,i2,i3) = LogHist(i1,i2,i3) + 1.
		enddo
	enddo

	do j = 1, NBin
		HistBinVal1(j) = EMin1 + EDelta1 * (real(j) - 0.5)
		HistBinVal2(j) = EMin2 + EDelta2 * (real(j) - 0.5)
		HistBinVal3(j) = EMin3 + EDelta3 * (real(j) - 0.5)
	enddo

	!take logs
	LogHist = log(LogHist)
	LogNFrame = log(real(NFrame))

	!create initial prob and free energy est
	FK = 0.

	!bin iteration to refine FK
	do i = 1, NIterBin
		do k = 1, NTemp
			do i1 = 1, NBin
				do i2 = 1, NBin
					do i3 = 1, NBin
						LogNum = LogHist(i1,i2,i3) - BetaK(k) * (WeightK1(k) * HistBinVal1(i1) &
							& + WeightK2(k) * HistBinVal2(i2) + WeightK3(k) * HistBinVal3(i3))
						LogDem = FK - BetaK * (WeightK1 * HistBinVal1(i1) &
							& + WeightK2 * HistBinVal2(i2) + WeightK3 * HistBinVal3(i3))
						LogDemMax = maxval(LogDem)
						LogP(i1,i2,i3) = LogNum - log(sum(exp(LogDem-LogDemMax))) - LogDemMax - LogNFrame
					enddo
				enddo
			enddo
			LogPMax = maxval(LogP)
			FK(k) = -log(sum(exp(LogP - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
	enddo
	!cleanup
	deallocate(LogHist)
	deallocate(LogP)
	deallocate(HistBinVal1, HistBinVal2)

	!iteration over all data (time consuming; ever so slightly more accurate)
	do i = 1, NIterAll
		do k = 1, NTemp
			do k2 = 1, NTemp
				do n = 1, NFrame
					E1 = EKN1(k2,n)
					E2 = EKN2(k2,n)
					E3 = EKN3(k2,n)
					LogNum = -BetaK(k) * (WeightK1(k) * E1 + WeightK2(k) * E2 + WeightK3(k) * E3)
					LogDem = FK - BetaK * (WeightK1 * E1 + WeightK2 * E2 + WeightK3 * E3)
					LogDemMax = maxval(LogDem)
					LogPKN(k2,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
				enddo
			enddo
			LogPMax = maxval(LogPKN)
			FK(k) = -log(sum(exp(LogPKN - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
	enddo

	!cleanup
	if (allocated(LogPKN)) deallocate(LogPKN)
end subroutine


subroutine log_weight_3(EKN1, EKN2, EKN3, WeightK1, WeightK2, WeightK3, &
	& BetaK, TargetWeight1, TargetWeight2, TargetWeight3, TargetBeta, &
	& FK, LogwKN, NTemp, NFrame)
	!This will calculate the log weights for general weights in the energy function.
	implicit none
	integer, intent(in) :: NTemp, NFrame
	real*8, dimension(NTemp, NFrame), intent(in) :: EKN1, EKN2, EKN3
	real*8, dimension(NTemp), intent(in) :: WeightK1, WeightK2, WeightK3
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, intent(in) :: TargetWeight1, TargetWeight2, TargetWeight3, TargetBeta
	real*8, dimension(NTemp), intent(in) :: FK
	real*8, dimension(NTemp, NFrame), intent(out) :: LogwKN
	real*8 :: LogNFrame
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax
	integer :: k, n
	!calculate individual probabilities
	LogNFrame = log(real(NFrame))
	do k = 1, NTemp
		do n = 1, NFrame
			LogNum = -TargetBeta * (EKN1(k,n) * TargetWeight1 + EKN2(k,n) * TargetWeight2 &
				& + EKN3(k,n) * TargetWeight3)
			LogDem = FK - BetaK * (EKN1(k,n) * WeightK1 + EKN2(k,n) * WeightK2 &
				& + EKN3(k,n) * WeightK3)
			LogDemMax = maxval(LogDem)
			LogwKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
		enddo
	enddo
	!shift probabilities
	LogwKN = LogwKN - maxval(LogwKN)	
end subroutine




!======== MULTIPLE TERMS IN THE ENERGY FUNCTION ========

subroutine free_energy_multi(EIKN, WeightIK, BetaK, MaxIter, RelTol, FK, NIter, NTerm, NTemp, NFrame)
	implicit none
	integer, intent(in) :: NTerm, NTemp, NFrame
	real*8, dimension(NTerm, NTemp, NFrame), intent(in) :: EIKN
	real*8, dimension(NTerm, NTemp), intent(in) :: WeightIK
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, intent(in) :: RelTol
	integer, intent(in) :: MaxIter
	real*8, dimension(NTemp), intent(out) :: FK
	integer, intent(out) :: NIter
	real*8, dimension(NTemp) :: OldFK
	real*8 :: LogNFrame
	real*8, dimension(:,:), allocatable :: LogPKN
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax, LogPMax, AbsTol
	integer :: i, j, k, n

	!make allocations
	allocate(LogPKN(NTemp, NFrame))

	!take logs
	LogNFrame = log(real(NFrame))

	!create initial prob and free energy est
	FK = 0.
	OldFK = 1.d200
	AbsTol = 0.

	!iteration over all data
	NIter = 0
	do while (NIter < MaxIter .and. any(abs(FK - OldFK) > AbsTol))
		OldFK = FK
		do i = 1, NTemp
			do k = 1, NTemp
				do n = 1, NFrame
					LogNum = -BetaK(i) * sum(WeightIK(:,i) * EIKN(:,k,n))
					do j = 1, NTemp
						LogDem(j) = FK(j) - BetaK(j) * sum(WeightIK(:,j) * EIKN(:,k,n))
					enddo 
					LogDemMax = maxval(LogDem)
					LogPKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
				enddo
			enddo
			LogPMax = maxval(LogPKN)
			FK(i) = -log(sum(exp(LogPKN - LogPMax))) - LogPMax
		enddo
		FK = FK - FK(1)
		!tolerance is relative to the maximum diff in FK
		AbsTol = RelTol * (maxval(FK) - minval(FK))
		NIter = NIter + 1
	enddo

	!cleanup
	if (allocated(LogPKN)) deallocate(LogPKN)
end subroutine


subroutine log_weight_multi(EIKN, WeightIK, BetaK, TargetWeightI, TargetBeta, FK, LogwKN, NTerm, NTemp, NFrame)
	!This will calculate the log weights for general weights in the energy function.
	implicit none
	integer, intent(in) :: NTerm, NTemp, NFrame
	real*8, dimension(NTerm, NTemp, NFrame), intent(in) :: EIKN
	real*8, dimension(NTerm, NTemp), intent(in) :: WeightIK
	real*8, dimension(NTemp), intent(in) :: BetaK
	real*8, dimension(NTerm), intent(in) :: TargetWeightI
	real*8, intent(in) :: TargetBeta
	real*8, dimension(NTemp), intent(in) :: FK
	real*8, dimension(NTemp, NFrame), intent(out) :: LogwKN
	real*8 :: LogNFrame
	real*8, dimension(NTemp) :: LogDem
	real*8 :: LogNum, LogDemMax
	integer :: k, n, j
	!calculate individual probabilities
	LogNFrame = log(real(NFrame))
	do k = 1, NTemp
		do n = 1, NFrame
			LogNum = -TargetBeta * sum(EIKN(:,k,n) * TargetWeightI)
			do j = 1, NTemp
				LogDem(j) = FK(j) - BetaK(j) * sum(EIKN(:,k,n) * WeightIK(:,j))
			enddo
			LogDemMax = maxval(LogDem)
			LogwKN(k,n) = LogNum - log(sum(exp(LogDem - LogDemMax))) - LogDemMax - LogNFrame
		enddo
	enddo
	!shift probabilities
	LogwKN = LogwKN - maxval(LogwKN)	
end subroutine