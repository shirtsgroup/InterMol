#!/usr/bin/python

VERSION=0
DATE="24 April 2007 1:25 pm"
AUTHOR="D. Ensign"

MODIFICATIONS_BY = "V.Voelz"
LAST_MODIFIED = "9/05/2007"


# HISTORY
# 9/05/07 VAV - Added plain Coulomb table (derived class of DielectricTable


# TO DO:
# * Implement more table types:
#
#  SF	- implement me!
#  SF2  - implement me!
#  DEL  - implement me!
#  DRF  - implement me!
#  SF3  - Shen Freed table 3 (?), implemented
#  Switch_Shift - GROMACS-style switch and shift for coulomb
#  Reaction_Field - RF with optional switch

#  Reaction_Field - RF with optional switch

usage_examples = """

# table=TableMaker(spacing=0.006) # otherwise use default parameters
# table=TableMaker(spacing=0.002)

# tabletype=DielectricTable()	
# tabletype=SF3(length_scale=4.0)
# tabletype=Switch_Shift()
# tabletype=ReactionField(rc=1.0)

# settings for eps->Infinity
# tabletype.krf = 1.0 / ( 2.0 * tabletype.rc )
# tabletype.crf = 3.0 / ( 2.0 * tabletype.rc ) 

# table=table.make(tabletype)

# print table
"""
	



class TableMaker:
	"""
	Class for making User-generated coulomb/vdw tables for use in GROMACS.

	Will require a table of some kind for the 'self.make' function; these tables are derived from 
	the 'DielectricTable' class below. 

	The function 'self.make' actually prints the table. The columns are:
		r	f	f''	g	g''	h	h''

		r   - the distance (ie, between two atoms for which we're computing the long-range 
			force)
		f   - the coulomb scaling rule
		f'' - second derivative of f
		g   - attractive (dispersive) part of LJ 
		g'' - second derivative of g
		h   - repulsive part of LJ
		h'' - second derivative of h

	From the manual, the long-range potential is 
	
			  qi qj 
		V(rij) = ------- f(rij)	+ C6 g(rij) + C12 h(rij)
			 4 pi e0

	with f(r), g(r), and h(r) defined above. Other values:
		
		rij - interatomic distance
		qi  - the charge of atom i
		pi  - pi, the value
		e0  - electrostatic permittivity, or permeability, whatsitcalled
		C6  - coefficient of r^6 term in the Lennard-Jones interaction
		C12 - coefficient of r^12 term in the Lennard-Jones interaction
	
	So you see we need not assume Lennard-Jones (6-12) at all. This is the default for the
	DielectricTable class, however.	 
 
	attributes:
	-----------
	spacing		distance between table points 		default 0.002 nm
	length		length of the table in nm 		default 5.0
	cutoff		distance at which to start entering	default 0.04
			nonzeros in the table

	methods:
	--------
	make		returns a string representing the table, ready to be written to file
	formatrow	used to print one row of the table 
	"""

	def __init__(self, spacing=0.002, length=50.0, cutoff=0.04):
		self.spacing=spacing
		self.length=length
		self.cutoff=cutoff
		
		# keep a list off all the derived classes of DielectricTable() we've added so far
		# This is useful for making decisions in System.prepare() 
		self.tabletype_names = []
		self.tabletype_names.append('Switch_Shift')
		self.tabletype_names.append('Coulomb')
		self.tabletype_names.append('ReactionField')
		self.tabletype_names.append('SF3')
			
		
	def make(self,tabletype):
		table=""
		r=0
		col1=col2=col3=col4=col5=col6=0
		while r<self.length:
			if r<self.cutoff:
				row=self.formatrow( r, col1, col2, col3, col4, col5, col6 )
			else:
				col1=tabletype.coulomb(r)
				col2=tabletype.coulomb_derivative2(r)
				col3=tabletype.dispersiveLJ(r)
				col4=tabletype.dispersiveLJ_derivative2(r)
				col5=tabletype.repulsiveLJ(r)
				col6=tabletype.repulsiveLJ_derivative2(r)	

				row=self.formatrow( r, col1, col2, col3, col4, col5, col6 )

			table += row
			r += self.spacing
		# the above has an extra /n at the end, so remove it
		table=table.rstrip()

		return table
	
	# this is ugly
	# could we consider just computing the *columns instead of the rows?*
	def formatrow( self, col0, col1, col2, col3, col4, col5 , col6, strfrm="% 10.5e  " ):
		row=(col0, col1, col2, col3, col4, col5 , col6 )
		rowstring=""
		for col in row:
			rowstring += strfrm % col
		rowstring += "\n"
		return rowstring
			
class DielectricTable:
	"""
	The base class for User-generated coulomb/vdw tables for use in GROMACS. 

	To make a new kind of table, make a new class and overwrite the functions here,
	if desired. IMPORTANT: double check with another method that you're getting the
	right numbers in the table. 
	
	methods:
	--------
        coulomb()			default return 1
        coulomb_derivative2()		default return 1
        dispersiveLJ()			default -1/pow(r,6)
        dispersiveLJ_derivative2()	default -42/pow(r,8) 
        repulsiveLJ()			default 1/pow(r,12)
        repulsiveLJ_derivative2()	default 156/pow(r,14)
	switchfunction()		not yet implemented

	attributes:
	-----------
	no general attributes
	"""
	
	def __init__( self ):
		pass

	def coulomb(self, r):
		return 1

	def coulomb_derivative2(self, r):
		return 1
	
	def dispersiveLJ(self, r):
		"""
		The traditional ~1/r^6 attractive LJ term.
		"""
		return -1/pow(r,6)
	
	def dispersiveLJ_derivative2(self, r):
		"""
		Second derivative of the traditional 1/r**6 attractive LJ Term.
		"""
		return -42/pow(r,8)
	
	def repulsiveLJ(self, r):
		"""
                The traditional ~1/r^12 repulsive LJ term.
                """
		return 1/pow(r,12)
	
	def repulsiveLJ_derivative2(self,r):
		"""
                Second derivative of the traditional 1/r**6 attractive LJ Term.
                """
		return 156/pow(r,14)

	def switchfunction( self ):
		pass

class Coulomb( DielectricTable ):
	""" GROMACS-style plain coulomb table for V = 1/r^a
	Here, a=1 is the standard coulomb potential, but we have generalized it to other exponents. 
	
	OPTIONAL INPUTS
	
	a        the 1/r^a exponent
        rc       the cutoff radius beyond which coulomb V(r) = 0 
	         (The default here is set to 100.0 nm, so there is effectively *no* cutoff

	"""

	def __init__( self, a=1, rc = 100. ):
		
		self.name = 'Coulomb'
		# Default is the Coulomb potential, with a=1
		self.a=float(a)
		self.rc=float(rc)       # the cutofff radius for the coulomb interaction
		                        # (default is LARGE, i.e, effectively no cutoff) 
		
		
	def coulomb( self, r ):
		if r<=self.rc:
			V=1/r**self.a
		else: V=0
		return V

	def coulomb_derivative2( self, r ):
		if r<=self.rc:
			D2V = (self.a+1.0)/r**(self.a+2.0) 
		else: D2V=0
	        return D2V




class Switch_Shift( DielectricTable ):
	""" GROMACS-style switch/shift

	Define three regions:
		<1>: 0 < r <= rs
		<2>: rs < r <= rc
		<3>: r > rc

	Define three constants:
		A = -( ( a+4 )*rc - ( a+1 )*rs ) / ( rc^(a+2) * ( rc-rs )^3 )
		B = ( ( a+3 )*rc - ( a+1 )*rs ) / ( rc^(a+2) * ( rc-rs )^3 )
		C = 1/rc^a - (A/3)*(rc-rs)^3 - (B/4)*(rc-rs)^4

	Where
		a = exponent for power law of the interaction, eg a=1 for coulomb: V(r) = 1/r
		rc = cutoff, after which potential/force/curvature = 0
		rs = switch radius, when to start applying the switch

	Then the potential, force, and curvature are given by:
		       { 1/r^a - C					, r in <1>
		V(r) = { 1/r^a - (A/3)*(r-rs)^3 - (B/4)*(r-rs)^4 - C	, r in <2>
		       { 0						, r in <3>

		       { 1/r^(a+1)				, r in <1> 
		F(r) = { 1/r^(a+1) + A(r-rs)^2 + B(r-rs)^3 	, r in <2>
		       { 0					, r in <3>

		       { (a+1)/r^(a+2)				, r in <1>
		R(r) = { (a+1)/r^(a+2) -2A(r-rs) -3B(r-rs)^2	, r in <2>
		       { 0					, r in <3>

	These functions are designed such that the potential resembles the unshifted potential at close range,
	but so that the potential, force, and curvature go smoothly to zero at rc. In addition, the potential
	and force are smooth at rs.

	attributes:
		a, default a=1
		rc, default rc=1.2
		rs, default rs=1.0
		A, B, and C as defined above

	WISH LIST: generalize the switching to work with any 1/r**a function, not just Coulomb. 

	"""

	def __init__( self, a=1, rc=1.2, rs=1.0 ):
		
		self.name = 'Switch_Shift'
		# Default is the Coulomb potential, with a=1
		self.a=float(a)
		self.rc=float(rc)
		self.rs=float(rs)
		self.A=self.calculateA()
		self.B=self.calculateB()
		self.C=self.calculateC()

	def calculateA( self ):
		a = self.a
		rc = self.rc
		rs = self.rs
		return -( ( a+4.0 )*rc - ( a+1.0 )*rs ) / ( rc**(a+2.0) * ( rc -rs )**3.0 )
		
	def calculateB( self ):
		a = self.a 
		rc = self.rc
		rs = self.rs
		return ( ( a+3.0 )*rc - ( a+1.0 )*rs ) / ( rc**(a+2.0) * ( rc-rs )**3.0 )

	def calculateC( self ):
		a = self.a
		rc = self.rc
		rs = self.rs
		A = self.A
		B = self.B
		return 1/rc**a - (A/3.0)*(rc-rs)**3.0 - (B/4.0)*(rc-rs)**4.0
		
	def coulomb( self, r ):
		if r<self.rs:
			V=1/r**self.a -self.C
		elif r<self.rc:
			V=1/r**self.a - self.C - (self.A/3.0 )*(r-self.rs)**3 - (self.B/4.0)*(r-self.rs)**4
		else: V=0
		return V

	def coulomb_derivative2( self, r ):
		if r<self.rs:
			D2V = (self.a+1.0)/r**(self.a+2.0) 
		elif r<self.rc:
			D2V = (self.a+1.0)/r**(self.a+2.0) - 2.0*self.A*(r-self.rs) - 3.0*self.B*(r-self.rs)**2
		else: D2V=0
		return D2V

class ReactionField( DielectricTable ):
	"""
	The reaction field equation is

		V(r) = (1/r)  + krf*r^2 -crf
	
	with:

		krf = (1/rc**3) * (eps - 1)/(2*eps + 1)
		crf = (1/rc) + krf*rc**2
		rc = cutoff, after which we assume a dielectric environment with eps, default 1.2 nm
		eps = the continuum dielectric constant, default = 78
	"""

	def __init__( self, eps=78, rc=1.2 ):
		self.name = 'ReactionField'
		self.eps=eps
		self.rc=rc
		self.krf=self.calculate_krf()
	 	self.crf=self.calculate_crf()

	def calculate_krf( self ):
		return ( 1.0/self.rc**3.0 ) * (self.eps - 1.0)/(2.0* self.eps + 1)

	def calculate_crf( self ):
		return ( 1.0/self.rc ) + self.krf * self.rc**2.0

	def coulomb( self, r ):
		if r < self.rc: V=self.krf * r**2 + 1/r - self.crf
		else: V=0
		#V=self.krf * r**2 + 1/r - self.crf
		return V
	
	def coulomb_derivative2( self, r ):
		if r < self.rc: Vpp = 2.0*self.krf + 2/r**3
		else: Vpp=0
		#Vpp = 2.0*self.krf + 2/r**3
		return Vpp

class SF3(DielectricTable):
	"""
	Shen-Freed dielectric table (number three?)

	functional form of dielectric constant:
	eps = Do + .5 (Do - Di) exp(-r s) ( (r s)^2 + 2 r s + 2 )

	methods:
	--------
	same as base class (DielectricTable)

        attributes:	formula symbol: description:				default values:
        -----------	--------------- ------------				---------------				
        epsilon_r       symbol Do	External dielectric constant            	78.0
        epsilon_inner   symbol Di	Inner (contact) dielectric constant      	 1.0
        length_scale    symbol s	characteristic length scale              	 3.0

	"""

	def __init__( self, epsilon_r=78.0, epsilon_inner=1.0, length_scale=3.0 ):

		self.name = 'SF3'		
		from math import exp as mathexp
                self.exp=mathexp

                self.epsilon_r=epsilon_r
                self.epsilon_inner=epsilon_inner
                self.length_scale=length_scale

                self.eps_avg= (epsilon_r - epsilon_inner) / 2.0
		DielectricTable.__init__(self)

        def coulomb(self,r):
		"""
		The Coulomb part is f(r) = 1 / (eps * r)
		"""
		eps=self.eps(r)
                return 1/(r * eps)

        def coulomb_derivative2(self, r):
		"""
		the derivative is
		 
		 d^2  [     1	 ]	2 f^2 + 2 r f f' + 2 r^2 (f')^2 - r^2 f f''
		 ---- [ -------- ]  =	-------------------------------------------
		 dr^2 [ r * f(r) ]			f^3 r^3
		
		It's fairly easy to do with the chain rule
		"""

		f=self.eps(r) 	# the dielectric funtion itself
		f2=pow(f,2)	# powers of the dielectric function
		f3=pow(f,3)

		r2=pow(r,2)	# powers of r
		r3=pow(r,3)

		fp = self.deps(r) # derivative of the dielectric function "f prime"
		fp2= pow(fp,2)    # derivative of the dielectric function, squared
		fpp= self.d2eps(r)# second derivative of the dielectric function, "f prime prime"

		num =  2 * f2  +  2 * f * r * fp  +  2 * r2 *fp2  -  f * r2 * fpp
		den =  f3 * r3		
		der2 = num / den
                return der2

	def eps(self, r):
		"""the dielectric function itself"""
	 	SR=self.length_scale*r
                eps=  self.epsilon_r - self.eps_avg * ( SR**2 + 2*SR  + 2 ) * self.exp(-SR)
		return eps

	def deps(self, r):
		"""the derivative of epsilon with respect to r"""
		SR=self.length_scale*r
		S3=pow(self.length_scale,3)
		R2=pow(r,2)
		deps = self.eps_avg * self.exp(-SR) * R2 * S3 
		return deps

	def d2eps(self, r):
		"""the second derivative of epsilon with respect to r"""
		SR=self.length_scale*r
		S3=pow(self.length_scale,3)
		d2eps= -self.eps_avg*self.exp(-SR)*r*S3*(-2+SR)
		return d2eps
	
