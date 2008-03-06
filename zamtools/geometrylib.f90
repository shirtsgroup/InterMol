
real*8 function ClipTrig(x)
	implicit none
	real*8, intent(in) :: x
	ClipTrig = min(1.D0, max(-1.D0, x))
end function

real*8 function NormRad(Rad)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: pi2 = 2.D0*pi
	real*8, intent(in) :: Rad
	NormRad = mod(Rad + pi, pi2) - pi
	if (NormRad < -pi) NormRad = NormRad + pi2
end function

real*8 function NormDeg(Deg)
	implicit none
	real*8, intent(in) :: Deg
	NormDeg = mod(Deg + 180.D0, 360.D0) - 180.D0
	if (NormDeg < -180.D0) NormDeg = NormDeg + 360.D0
end function


subroutine Centroid(Pos, Ret, N, Dim)
	implicit none
	integer, intent(in) :: N, Dim
	real*8, dimension(N,Dim), intent(in) :: Pos
	real*8, dimension(Dim), intent(out) :: Ret
	Ret = sum(Pos, 1) / real(N)
end subroutine

subroutine CentroidInd(Pos, AtomInd, Ret, N, Dim, NAtomInd)
	implicit none
	integer, intent(in) :: N, Dim, NAtomInd
	real*8, dimension(N,Dim), intent(in) :: Pos
	integer, dimension(NAtomInd), intent(in) :: AtomInd
	real*8, dimension(Dim), intent(out) :: Ret
	integer :: i
	Ret = 0.
	do i = 1, NAtomInd
		Ret = Ret + Pos(AtomInd(i) + 1, :)
	enddo
	Ret = Ret / real(NAtomInd)
end subroutine

real*8 function Length(Vec, Dim)
	implicit none
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: Vec
	Length = sqrt(sum(Vec*Vec))
end function

subroutine UnitVec(Vec, Ret, Dim)
	implicit none
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: Vec
	real*8, dimension(Dim), intent(out) :: Ret
	Ret = Vec / sqrt(sum(Vec*Vec))
end subroutine

subroutine CrossProd3(r1, r2, Ret)
	real*8, dimension(3), intent(out) :: Ret
	real*8, dimension(3), intent(in) :: r1, r2
	Ret = (/r1(2)*r2(3)-r1(3)*r2(2), r1(3)*r2(1)-r1(1)*r2(3), r1(1)*r2(2)-r1(2)*r2(1)/)
end subroutine

real*8 function RadFromTrig(SinVal, CosVal)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, intent(in) :: SinVal, CosVal
	real*8, external :: ClipTrig
	RadFromTrig = acos(ClipTrig(CosVal))
	if (SinVal < 0.) RadFromTrig = 2.*pi - RadFromTrig
end function

real*8 function DegFromTrig(SinVal, CosVal)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	real*8, intent(in) :: SinVal, CosVal
	real*8, external :: ClipTrig, RadFromTrig
	DegFromTrig = DegPerRad * RadFromTrig(SinVal, CosVal)
end function

real*8 function Angle3Rad(Pos1, Pos2, Pos3)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	real*8, dimension(3), intent(in) :: Pos1, Pos2, Pos3
	real*8, dimension(3) :: Vec21, Vec23
	real*8 :: Norm, Phi
	real*8, external :: ClipTrig, NormRad
	if (all(Pos1 == Pos2) .or. all(Pos2 == Pos3)) then
		Angle3Rad = 0.
		return
	endif
	Vec21 = Pos1 - Pos2
	Vec23 = Pos3 - Pos2
	Norm = sqrt(sum(Vec21*Vec21)*sum(Vec23*Vec23))
	Phi = ClipTrig(dot_product(Vec21, Vec23) / Norm)
	Phi = acos(Phi)
	Angle3Rad = NormRad(Phi)
end function

real*8 function Angle3Deg(Pos1, Pos2, Pos3)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	real*8, dimension(3), intent(in) :: Pos1, Pos2, Pos3
	real*8, external :: Angle3Rad
	Angle3Deg = Angle3Rad(Pos1, Pos2, Pos3) * DegPerRad
end function

real*8 function Dihedral3Rad(Pos1, Pos2, Pos3, Pos4)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	real*8, dimension(3), intent(in) :: Pos1, Pos2, Pos3, Pos4
	real*8, dimension(3) :: Vec12, Vec23, Vec34, Norm12, Norm34
	real*8 :: Norm, Phi
	real*8, external :: ClipTrig, NormRad
	external :: CrossProd3
	Vec12 = Pos2 - Pos1
	Vec23 = Pos3 - Pos2
	Vec34 = Pos4 - Pos3
	call CrossProd3(Vec12, Vec23, Norm12)
	call CrossProd3(Vec23, Vec34, Norm34)
	Norm = sqrt( sum(Norm12*Norm12) * sum(Norm34*Norm34) )
	Phi = ClipTrig(dot_product(Norm12, Norm34) / Norm)
	Phi = acos(Phi)
	if (dot_product(Vec12, Norm34) < 0.) Phi = -Phi
	Dihedral3Rad = NormRad(Phi)
end function

real*8 function Dihedral3Deg(Pos1, Pos2, Pos3, Pos4)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	real*8, dimension(3), intent(in) :: Pos1, Pos2, Pos3, Pos4
	real*8, external :: Dihedral3Rad
	Dihedral3Deg = Dihedral3Rad(Pos1, Pos2, Pos3, Pos4) * DegPerRad
end function

subroutine GetVecMappingRad(Vec1, Vec2, Vec, Ang, Dim)
	implicit none
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: Vec1, Vec2
	real*8, dimension(Dim), intent(out) :: Vec
	real*8, intent(out) :: Ang
	real*8, dimension(Dim) :: v1, v2
	real*8 :: CosAng
	real*8, external :: ClipTrig
	external :: CrossProd3, UnitVec
	call UnitVec(Vec1, v1, Dim)
	call UnitVec(Vec2, v2, Dim)
	CosAng = ClipTrig(dot_product(v1, v2))
	Ang = acos(CosAng)
	if (CosAng == 1.) then
		Vec = v1
	elseif (CosAng == -1.) then
		Vec(2:Dim) = v1(1:Dim-1)
		Vec(1) = v1(Dim)
		Vec = Vec - dot_product(Vec, v1)
	else
		call CrossProd3(v1, v2, Vec)
	endif
end subroutine

subroutine GetVecMappingDeg(Vec1, Vec2, Vec, Ang, Dim)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: Vec1, Vec2
	real*8, dimension(Dim), intent(out) :: Vec
	real*8, intent(out) :: Ang
	external :: GetVecMappingRad
	call GetVecMappingRad(Vec1, Vec2, Vec, Ang, Dim)
	Ang = Ang * DegPerRad
end subroutine

subroutine RotateArrayAboutPoint(Pos, RotMat, Point, Ret, N, Dim)
	integer, intent(in) :: N, Dim
	real*8, dimension(N, Dim), intent(in) :: Pos
	real*8, dimension(Dim, Dim), intent(in) :: RotMat
	real*8, dimension(Dim), intent(in) :: Point
	real*8, dimension(N, Dim), intent(out) :: Ret
	integer :: i
	do i = 1, N
		Ret(i,:) = Pos(i,:) - Point
	enddo
	Ret = matmul(Ret, RotMat)
	do i = 1, N
		Ret(i,:) = Ret(i,:) + Point
	enddo
end subroutine

subroutine RotatePointAboutPoint(Pos, RotMat, Point, Ret, Dim)
	integer, intent(in) :: Dim
	real*8, dimension(Dim), intent(in) :: Pos
	real*8, dimension(Dim, Dim), intent(in) :: RotMat
	real*8, dimension(Dim), intent(in) :: Point
	real*8, dimension(Dim), intent(out) :: Ret
	Ret = matmul(Pos - Point, RotMat) + Point
end subroutine


subroutine RotMat3Rad(Vec, Ang, Ret)
	implicit none
	real*8, dimension(3), intent(in) :: Vec
	real*8, intent(in) :: Ang
	real*8, dimension(3,3), intent(out) :: Ret
	real*8, dimension(3) :: v
	real*8 :: rcos, rsin
	external :: UnitVec
	Ret = 0.
	rcos = cos(Ang)
	rsin = sin(Ang)
	call UnitVec(Vec, v, 3)
	Ret(1,1) =         rcos + v(1)*v(1)*(1-rcos)
	Ret(2,1) =  v(3) * rsin + v(2)*v(1)*(1-rcos)
	Ret(3,1) = -v(2) * rsin + v(3)*v(1)*(1-rcos)
	Ret(1,2) = -v(3) * rsin + v(1)*v(2)*(1-rcos)
	Ret(2,2) =         rcos + v(2)*v(2)*(1-rcos)
	Ret(3,2) =  v(1) * rsin + v(3)*v(2)*(1-rcos)
	Ret(1,3) =  v(2) * rsin + v(1)*v(3)*(1-rcos)
	Ret(2,3) = -v(1) * rsin + v(2)*v(3)*(1-rcos)
	Ret(3,3) =         rcos + v(3)*v(3)*(1-rcos)
end subroutine

subroutine RotMat3Deg(Vec, Ang, Ret)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: RadPerDeg = pi / 180.D0
	real*8, dimension(3), intent(in) :: Vec
	real*8, intent(in) :: Ang
	real*8, dimension(3,3), intent(out) :: Ret
	external :: RotMat3Rad
	call RotMat3Rad(Vec, Ang * RadPerDeg, Ret)
end subroutine


subroutine RotMat3EulerRad(Phi, Theta, Psi, Ret)
	implicit none
	real*8, dimension(3,3), intent(out) :: Ret
	real*8, intent(in) :: Phi, Theta, Psi
	real*8 :: Sin1, Cos1, Sin2, Cos2, Sin3, Cos3
	Sin1 = sin(Phi)
	Cos1 = cos(Phi)
	Sin2 = sin(Theta)
	Cos2 = cos(Theta)
	Sin3 = sin(Psi)
	Cos3 = cos(Psi)
	Ret(1,1) = Cos1*Cos3 - Sin1*Cos2*Sin3
	Ret(1,2) = Sin1*Cos3 + Cos1*Cos2*Sin3
	Ret(1,3) = Sin2*Sin3
	Ret(2,1) = -Cos1*Sin3 - Sin1*Cos2*Cos3
	Ret(2,2) = -Sin1*Sin3 + Cos1*Cos2*Cos3
	Ret(2,3) = Sin2*Cos3
	Ret(3,1) = Sin1*Sin2
	Ret(3,2) = -Cos1*Sin2
	Ret(3,3) = Cos2
end subroutine

subroutine RotMat3EulerDeg(Phi, Theta, Psi, Ret)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: RadPerDeg = pi / 180.D0
	real*8, dimension(3,3), intent(out) :: Ret
	real*8, intent(in) :: Phi, Theta, Psi
	external :: RotMat3EulerRad
	call RotMat3EulerRad(Phi * RadPerDeg, Theta * RadPerDeg, Psi * RadPerDeg, Ret)
end subroutine


subroutine RotMat3Q(q0, q1, q2, q3, Ret)
	implicit none
	real*8, dimension(3,3), intent(out) :: Ret
	real*8, intent(in) :: q0, q1, q2, q3
	Ret(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
	Ret(1,2) = 2.*(q1*q2 + q0*q3)
	Ret(1,3) = 2.*(q1*q3 - q0*q2)
	Ret(2,1) = 2.*(q1*q2 - q0*q3)
	Ret(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
	Ret(2,3) = 2.*(q2*q3 + q0*q1)
	Ret(3,1) = 2.*(q1*q3 + q0*q2)
	Ret(3,2) = 2.*(q2*q3 - q0*q1)
	Ret(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3
end subroutine


subroutine EulerFromRotMat3Rad(RotMat, Phi, Theta, Psi)
	implicit none
	real*8, dimension(3,3), intent(in) :: RotMat
	real*8, intent(out) :: Phi, Theta, Psi
	real*8 :: sphi
	real*8, external :: RadFromTrig
	Theta = acos(min(1., max(-1., RotMat(3,3))))
	if (RotMat(3,3)==1. .or. (RotMat(3,1)==0. .and. RotMat(3,2)==0.) &
		& .or. (RotMat(1,3)==0. .and. RotMat(2,3)==0.)) then
		Psi = 0.
		Phi = RadFromTrig(RotMat(1,2), RotMat(1,1))
	else
		sphi = sqrt(1.-RotMat(3,3)*RotMat(3,3))
		Phi = RadFromTrig(RotMat(3,1)/sphi, -RotMat(3,2)/sphi)
		Psi = RadFromTrig(RotMat(1,3)/sphi, RotMat(2,3)/sphi)
	endif
end subroutine

subroutine EulerFromRotMat3Deg(RotMat, Phi, Theta, Psi)
	implicit none
	real*8, parameter :: pi = 3.1415926535897931D0
	real*8, parameter :: DegPerRad = 180.D0/pi
	real*8, dimension(3,3), intent(in) :: RotMat
	real*8, intent(out) :: Phi, Theta, Psi
	external :: EulerFromRotMat3Rad
	call EulerFromRotMat3Rad(RotMat, Phi, Theta, Psi)
	Phi = Phi * DegPerRad
	Theta = Theta * DegPerRad
	Psi = Psi * DegPerRad
end subroutine

