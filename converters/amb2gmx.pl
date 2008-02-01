#!/usr/bin/perl -w

# amb2gmx.pl v.001 2006-01-14
# ==============================================================================================
# A script to convert AMBER prmtop and coordinate files to gromacs topology and coordinate files.
# Originally written by Eric Sorin, Pande lab, Stanford.
# Modified by David Mobley and John Chodera, Dill lab, UCSF.
# ==============================================================================================
# Insert license here.
# ==============================================================================================
# REQUIREMENTS:
# - From AMBER:
#   - rdparm
# ==============================================================================================
# CHANGELOG:
# v.001 2006-01-14
# - Changelog initiated.
# v.002 2006-06-08 DLM
# - Modifying Qmax routine to enforce integer charge (nearest integer) rather than charge 
#   neutrality. This is done in the same way as previously -- modifying the atom with largest charge.
# v.003 2007-04-03 DLM
# - Modifying bond parsing routine for rdparm output; it used incorrect column widths and would fail when
#   bond constants got larger than kB=999.99. 
# v.004 2007-09-03 JDC
# - Added elementary support for truncated octahedron boxes. An error is thrown if box is not
#   rectilinear or truncated octahedron.
# v.005 2008-01-31 JDC
# - Hacked in support for magnesium ions.
# - Added support for writing .g96 files.
# - Fixed bug in truncated octahedron box support.  Things now compare favorably with AMBER energies.
# ==============================================================================================
# TODO and KNOWN BUGS:
# - Figure out what the "[ pairs ]" section is really doing.
# ==============================================================================================

#===============================================================================================
# REFERENCES:
#===============================================================================================
# If you use this in a publication, please reference: D. L. Mobley, J. D. Chodera, K. A. Dill. "On the use of orientational restraints and symmetry number corrections in alchemical free energy calculations", Journal of Chemical Physics 125:084902, 2006

# ==============================================================================================
# USE declarations
# ==============================================================================================

use strict;                 # Always, ALWAYS use strict!
use IO::File;               # For creating file handles in variables.
use Getopt::Long;           # For command-line argument processing.

#==============================================================================
# CONSTANTS
#==============================================================================

# converting ligand rdparm/pdb files to gro/top files
# dihedrals converted to RB functions...
my $cal = 4.184; # J/cal

#==============================================================================
# SUBROUTINES
#==============================================================================

# Trim leading and training spaces from a string, returning the result.
# @param arg the string to be trimmed
# @returns the string, with leading and trailing whitespace removed
sub trim {
  # get arguments 
  my $arg = shift @_;
  
  # remove leading and trailing spaces
  $arg =~ s/^\s+//;
  $arg =~ s/\s+$//;
  
  # return result
  return $arg;
}

# Determine whether a residue name should be considered 'solute' or 'solvent'.
# @param resname the residue name
# @returns 1 if solute, 0 if solvent.
sub is_solute {
  # get arguments
  my $resname = shift @_;
  
  if($resname eq "WAT" || $resname eq "Na+" || $resname eq "Cl-" || $resname eq "MG2") {
    return 0;
  }
  
  return 1;
}

#==============================================================================
# MAIN
#==============================================================================

# Parse command-line arguments.
my $prmtop_filename = "";
my $crd_filename = "";
my $outname = "";                 # prefix to be used for all rdparm and gromacs files
my $help_flag = 0;                # 1 if help has been requested
my $debug_flag = 0;               # 1 if debug information is to be printed
GetOptions( "prmtop=s" => \$prmtop_filename, "crd=s" => \$crd_filename, "outname=s" => \$outname, "debug" => \$debug_flag, "help" => \$help_flag );

# Make sure required arguments are specified.
if( ($prmtop_filename eq "") || ($crd_filename eq "") || ($outname eq "") ) {
  # Show help instead.
  $help_flag = 1;
}

# Display help if requested.
if( $help_flag ) {
  print "amb2gmx.pl\n";
  print "Convert an AMBER system to gromacs format.\n";
  print "\n";
  print "usage: ./amb2gmx.pl --prmtop amber-prmtop-filename --crd amber-crd-filename --outname outname [--debug] [--help]\n";
  print "\n";
  print "example: ./amb2gmx.pl --prmtop amber/1LMB_25wat96.top --crd amber/inpcrd --outname lambda\n";
  print "\n";
  exit(0);
}

# Make sure other required commands are present.
# TODO: Test for existence of rdparm in path.

# Create specially-formatted PDB file.
my $pdb_filename = "${outname}.pdb";
{
  my $command = "cat $crd_filename | ambpdb -atm -p $prmtop_filename > $pdb_filename";
  print "Calling rdparm to extract information from prmtop file.\n$command\n" if($debug_flag);
  my $output = `$command`;
  print "$output\n" if($debug_flag);
}

# Create input file for rdparm.
my $rdparm_input_filename = "${outname}.rdparm.in";
{
  print "Creating rdparm input file ${rdparm_input_filename}...\n" if($debug_flag);
  my $rdparm_input_file = new IO::File;
  open $rdparm_input_file, ">", $rdparm_input_filename;
  print $rdparm_input_file <<EOF;
atoms
bonds
angles
dihedrals
pertbonds
pertangles
pertdihedrals
printExcluded
printLennardJones
printTypes
info
quit
EOF
  close $rdparm_input_file;
}

# Use rdparm to extract required information from prmtop file.
my $rdparm_output_filename = "${outname}.rdparm.out";
{
  my $command = "rdparm $prmtop_filename < $rdparm_input_filename > $rdparm_output_filename";
  print "Calling rdparm to extract information from prmtop file.\n$command\n" if($debug_flag);
  my $output = `$command`;
  print "$output\n" if($debug_flag);
}

################  general I/O  ####################

open(PDB,$pdb_filename) || die "Error: can't open PDB file\n\n";
open(TXT,$rdparm_output_filename) || die "Error: can't open RDPARM file\n\n";
my $grofile = "$outname".".gro";
open(GRO,">$grofile") || die "Error: can't write to .gro file\n\n";
my $g96file = "$outname".".g96";
open(G96,">$g96file") || die "Error: can't write to .g96 file\n\n";
my $topfile = "$outname".".top";
open(TOP,">$topfile") || die "Error: can't write to .top file\n\n";

# Extract box information from AMBER CRD file.
# TODO: This relies on the coordinate file missing the velocity information -- fix this.
my $lastline = `tail -n 1 $crd_filename`;
# Store box extents, in nm.
my %box = ();
$box{x} = substr($lastline,0,12) * 0.1;
$box{y} = substr($lastline,12,12) * 0.1;
$box{z} = substr($lastline,24,12) * 0.1;
# store box angles, in degrees.
my %box_angles = ();
$box_angles{1} = substr($lastline,36,12);
$box_angles{2} = substr($lastline,48,12);
$box_angles{3} = substr($lastline,60,12);

#############	parse the pdb file  ###############	
print "Parsing PDB file $pdb_filename...\n";
my $pdbnumatoms = 0; # number of atoms in the pdb file
my %coordinates = (); # hash containing coordinates.  "${atomindex}x" is x-coordinates of atom $atomindex, which runs from 1...natoms.
my %resname = (); # hash containing residue names corresponding to atom numbers.
my %resnum = (); # hash contining residue numbers corresponding to atom numbers.
my %pdbname = (); # hash containing name of atom in PDB file

my $firstsodium = -1; # residue number of first sodium ion
my $firstchloride = -1; # residue number of first chloride ion
my $firstmagnesium = -1; # residue number of first magnesium
my $firstwat = -1; # residue number of first water

my $nproteinres = 0;
my $nsodium = 0;
my $nchloride = 0;
my $nmagnesium = 0;
my $nwat = 0;

my $last_resnum = -1;

while(my $line = <PDB>) {
  chomp $line;

  # Skip REMARK and TER fields.
  next if(substr($line,0,6) eq "REMARK" || substr($line,0,6) eq "TER   " || substr($line,0,6) eq "END   ");
  
  # Parse line according to specification from ambpdb for -atm flag:
  #    50 FORMAT(3f10.3,3i5,5x,a4,i5,1x,a4)
  #            write(nf,50) (c(jc+m),m=1,3),
  #   .          itype,isrn,nn,lbres(j),j+ioffset,igraph(k)

  # Extract x,y,z and convert to nm.
  my $x = substr($line,0,10)  * 0.1;
  my $y = substr($line,10,10) * 0.1;
  my $z = substr($line,20,10) * 0.1;
  
  # Extract atom type index.
  my $itype = substr($line,30,5) + 0;

  # TODO: I am not sure what 'isrn' and 'nn' represent -- find out!
  my $isrn = substr($line,35,5) + 0;
  my $nn = substr($line,40,5) + 0;

  # Extract residue name.
  my $this_resname = trim(substr($line,50,4));

  # Extract residue number.
  my $this_resnum = substr($line,54,5) + 0;

  # Extract atom name.
  my $atomname = substr($line,60,4);
  
  # print "x = $x\ny = $y\nz = $z\nitype = $itype\nisrn = $isrn\nnn = $nn\nlbres = $lbres\nindex = $index\nigraph = $igraph\n";

  # Store.
  $pdbnumatoms++;
  $resname{$pdbnumatoms} = $this_resname;
  $resnum{$pdbnumatoms} = $this_resnum;
  $pdbname{$pdbnumatoms} = $atomname;
  $coordinates{"${pdbnumatoms}x"} = $x;
  $coordinates{"${pdbnumatoms}y"} = $y;
  $coordinates{"${pdbnumatoms}z"} = $z;

  # Make a note of the first sodium, chloride, magnesium, and solvent residues.
  # TODO: Generalize this to track multiple types of arbitrary molecules?
  if($firstsodium == -1 && $this_resname eq "Na+") { $firstsodium = $this_resnum; }
  if($firstchloride == -1 && $this_resname eq "Cl-") { $firstchloride = $this_resnum; }
  if($firstmagnesium == -1 && $this_resname eq "MG2") { $firstmagnesium = $this_resnum; }
  if($firstwat == -1 && $this_resname eq "WAT") { $firstwat = $this_resnum; }

  # accumulate counts if this is a new residue
  if($this_resnum != $last_resnum) {
    if($this_resname eq "WAT") {
        $nwat++;
    } elsif($this_resname eq "Na+") {
        $nsodium++;
    } elsif($this_resname eq "Cl-") {
        $nchloride++;
    } elsif($this_resname eq "MG2") {
        $nmagnesium++;
    } else {
        $nproteinres++;
    }
  }

  $last_resnum = $this_resnum;
}
close(PDB);

# Total number of residues.
my $ntotalres = $resnum{$pdbnumatoms};

printf STDOUT "found:\n";
printf STDOUT "%5d protein or solute residues\n", $nproteinres;
printf STDOUT "%5d sodium ions\n", $nsodium;
printf STDOUT "%5d chloride ions\n", $nchloride;
printf STDOUT "%5d magnesium ions\n", $nmagnesium;
printf STDOUT "%5d waters\n", $nwat;
printf STDOUT "\n";
printf STDOUT "%5d TOTAL residues read\n", $ntotalres;

#############	read coordinates from AMBER .crd file	###############	

{
    print "Reading coordinates from CRD file...";
    open(CRD,$crd_filename) || die "Error: can't open CRD file\n\n";
    # title line
    my $line = <CRD>;
    # atom number line
    $line = <CRD>;
    # read atom numbers
    my @coordinate_stack = ();
    while($#coordinate_stack+1 < 3*$pdbnumatoms) {
	$line = <CRD>;
	chomp $line;
	my @elements = split " ", $line;
	push @coordinate_stack, @elements;
    }
    close(CRD);
    for(my $atomindex = 1; $atomindex <= $pdbnumatoms; $atomindex++) {
	$coordinates{"${atomindex}x"} = 0.1 * shift @coordinate_stack;
	$coordinates{"${atomindex}y"} = 0.1 * shift @coordinate_stack;
	$coordinates{"${atomindex}z"} = 0.1 * shift @coordinate_stack;
    }
}

#############	parse the prmtop file	###############	
print "Reading prmtop file...\n";
# TODO: This only works with AMBER7 format top files.
my %atom_type_indices = ();
my $atomindex = 0;
open(PRMTOP,$prmtop_filename) || die "Error: can't open prmtop file \"$prmtop_filename\"\n\n";
while(my $line = <PRMTOP>) {
    chomp $line;
    
    if(substr($line,0,21) eq "%FLAG ATOM_TYPE_INDEX") {
	# eat format line
	<PRMTOP>;
	# read atom type indices
	while($atomindex < $pdbnumatoms) {
	    # read a line
	    $line = <PRMTOP>;	    
	    chomp $line;
	    my $column_index = 0;
	    while($column_index < length $line) {
		$atomindex++;
		# parse I8 field
		my $type_index = substr($line,$column_index,8) + 0;
		$atom_type_indices{$atomindex} = $type_index;
		$column_index += 8;
	    }
	}
    }
}
close(PRMTOP);
print "$atomindex atom type indices read.\n\n";

#############	parse the dump file	###############	
my $nato = 0; my $getatoms=0;
my $nbon = 0; my $getbonds=0;
my $npai = 0; my $getpairs=0;
my $nang = 0; my $getangs =0;
my $ndih = 0; my $getdihs =0;
my $nimp = 0; my $getimps =0;
my $nlj  = 0; my $getlj   =0;

# TODO: Make this a structure?

# atom type information
my %gaff;     # atom type for each atom
my %atomname; # atom name
my %charge;   # charge
my %mass;     # mass (in amu)

# bond type information
my %bstr;     # spring constant
my %leng;     # equilibrium length
my (%atom1, %atom2); # atom indices

# angle information
my %str;      # angle force constant, in cal / ???
my %ang;      # equilibrium angle
my (%atmo1, %atmo2, %atmo3);  # atom indices

# torsion information
my %idivf;
my %pk;
my %phase;
my %pn;
my %at;

# parilist information
my (%pair1, %pair2);

# LJ type information
my %rad; # LJ radius (sigma)
my %eps; # LJ well depth (epsilon)

while(my $line = <TXT>) {
  chomp $line;

  my $splitline = $line;
  for($splitline) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
  my @input = split(/ /,$splitline);

  # skip empty lines
  next if( length($line)==0 || $#input<0 );

  # get atom names, charges, masses and atom types
  if(($input[0] eq 'RDPARM')&&($getatoms==1)) { 
    $getatoms=0; 
    next;
  } elsif($getatoms==1) {
    # parse line if it contains an atom entry
    if(substr($line,0,6) =~ /\d+/) {
      # parse line
      my $this_atomnumber = substr($line,0,6) + 0;
      my $this_atomname   = trim(substr($line,9,4));
      my $this_charge     = substr($line,14,8) + 0.0;
      my $this_mass       = substr($line,23,4) + 0.0;
      my $this_resnum     = substr($line,29,4) + 0;
      my $this_resname    = trim(substr($line,34,4));
      my $this_atomtype   = trim(substr($line,41,3));
      my $this_tree       = trim(substr($line,47,3));

      # store
      $nato++;
      $atomname{$nato} = $this_atomname;
      $charge{$nato}   = $this_charge;
      $mass{$nato}     = $this_mass;
      $gaff{$nato}     = $this_atomtype;

      #if($nato < 100) { print STDOUT "$nato $atomname{$nato} $charge{$nato} $mass{$nato} $gaff{$nato}\n"; }
    }
  } elsif($input[0] eq 'Number:') { 
    $getatoms = 1; 
  }
   
  # get bonds #
  if(($input[0] eq 'RDPARM')&&($getbonds==1)) {
    $getbonds=0; 
    next;
  } elsif($getbonds==1){  

    # parse line if it contains a bond entry
    if(substr($line,0,8) =~ /\d+/) {
      # parse line
      my $this_bondindex = substr($line,0,8) + 0;
      my $this_kB = substr($line,9,9) * 2 * 100 * $cal; # in kJ/mol/nm^2?
      my $this_Req = substr($line,19,6) * 0.1; # in nm
      my $rest_of_line = substr($line,25);
      for($rest_of_line) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
      my ($this_atomname1, $this_atomname2, $atomnumbers) = split " ", $rest_of_line;
      # extract atomnumbers from (x,y) format.
      $atomnumbers =~ /\((\d+),(\d+)\)/;
      my $this_atomnumber1 = $1; my $this_atomnumber2 = $2;
      
      # store
      $nbon++;
      $bstr{$nbon} = $this_kB;
      $leng{$nbon} = $this_Req;
      $atom1{$nbon} = $this_atomnumber1;
      $atom2{$nbon} = $this_atomnumber2;

      # if($nbon < 100) { print STDOUT "$nbon $bstr{$nbon} $leng{$nbon} $atom1{$nbon} $atom2{$nbon}\n"; }
    }
  } elsif($input[0] eq 'Bond'){ 
    $getbonds = 1; 
  }

  # get angles #
  if(($input[0] eq 'RDPARM')&&($getangs==1)){ $getangs=0; next;
  }elsif($getangs==1){  
    my $test = $input[0];
    for($test){ s/\://; }
    if($test > 0){
      $nang++;
      $str{$nang}  = $input[1]*2*$cal;
      $ang{$nang}  = $input[2];
      my @bonded = split(/\,/,$input[6]);
      $atmo1{$nang} = $bonded[0];
      for($atmo1{$nang}){ s/\(//; }
      $atmo2{$nang} = $bonded[1];
      $atmo3{$nang} = $bonded[2];
      for($atmo3{$nang}){ s/\)//; }
      # print STDOUT "$nang $str{$nang} $ang{$nang} $atmo1{$nang} $atmo2{$nang} $atmo3{$nang}\n";
    }
  }elsif($input[0] eq 'Angle'){ $getangs = 1; }


  # get torsions and take pairs from 1-4 interactions #
  if(($input[0] eq 'RDPARM')&&($getdihs==1)){ $getdihs=0; next;
  }elsif($getdihs==1){  
    my $test = $input[0];
    for($test){ s/\://; }
    my $addtopairlist=1;  
    if(($test eq "E")||($test eq "B")){ 
      for(my $in=0;$in<=8;$in++){
        $input[$in]=$input[$in+1]
      }
      $test = $input[0];
      for($test){ s/\://; }
      $addtopairlist=0;
    }

    if($test > 0){
      $ndih++;
      $idivf{$ndih} = 1; 
      $pk{$ndih}    = $input[1];
      $phase{$ndih} = $input[2];
      $pn{$ndih}    = $input[3];

      # JDC: This appears to do nothing, so I removed it.
      #$quartet = $input[8];
      #if($quartet eq $quartet_old){ }
      #$quartet_old = $quartet;

      my @bonded = split(/\,/,$input[8]);
      $at{$ndih,1} = $bonded[0];
      for($at{$ndih,1}){ s/\(//; }
      $at{$ndih,2} = $bonded[1];
      $at{$ndih,3} = $bonded[2];
      $at{$ndih,4} = $bonded[3];
      for($at{$ndih,4}){ s/\)//; }
      if ($addtopairlist==1) {                  
              $npai++;                          
              $pair1{$npai}= $at{$ndih,1};      
              $pair2{$npai}= $at{$ndih,4};      
      }                                        
    }
  }elsif($input[0] eq 'Dihedral'){ $getdihs = 1; }

  # get LJ params #
  if((substr($line,0,6) eq 'RDPARM')&&($getlj==1)){ $getlj=0; next;
  }elsif($getlj==1){
    # parse line if it contains a bond entry
    if(substr($line,9,6) =~ /\d+/) {
      # parse line
      my $atomtype = trim(substr($line,3,3));
      my $radius   = substr($line,9,6) * (2.**(5./6.))/10.;
      my $epsilon  = substr($line,17,6) * $cal;

      # store
      $nlj++;
      $rad{$nlj} = $radius;
      $eps{$nlj} = $epsilon;

      printf STDOUT "LJ %3d %6.4f %6.4f\n", $nlj, $radius, $epsilon;
    }
  }elsif($input[0] eq 'Type'){ $getlj = 1; }
}

print STDOUT "found:
	$nato atoms
	$nbon bonds
	$nang angles
	$ndih torsions\n";
	
my $numpropers = $ndih;
close(TXT);


################  print gro/top headers  ####################
my $date = `date`; chomp $date;

print TOP "; $outname.top created by rdparm2gmx.pl $date

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333\n";

print GRO "$outname.gro created by rdparm2gmx.pl $date\n $nato\n";


#############	write .gro file		 ###############	
if($pdbnumatoms != $nato){ 
  print STDOUT "\n\tError: pdb and rdparm disagree on number of atoms!\n\n";
  exit();
}else{
  for(my $i=1;$i<=$nato;$i++) {
        printf GRO "%5d%-4s%6s%5d%8.3f%8.3f%8.3f\n",$resnum{$i},$resname{$i},$atomname{$i},$i,$coordinates{"${i}x"},$coordinates{"${i}y"},$coordinates{"${i}z"};

  }
  # print GRO "   10.00000   10.00000   10.00000\n";
  # rectangular box
  if($box_angles{1} == "90.0000000" && $box_angles{2} == "90.0000000" && $box_angles{3} == "90.0000000"){
    print "Rectangular box.";
    printf GRO "%11.5f%11.5f%11.5f\n", $box{x}, $box{y}, $box{z};
  } elsif($box_angles{1} == "109.4712190" && $box_angles{2} == "109.4712190" && $box_angles{3} == "109.4712190") {
    print "Octahedral box.";
    # TODO: Check to make sure all box dimensions are equal.
    my $d = $box{x} + 0.0;
    printf GRO "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n", $d, (2./3.)*sqrt(2)*$d, (1./3.)*sqrt(6)*$d, 0.0, 0.0, (1.0/3.0)*$d, 0.0, -(1./3.)*$d, (1./3.)*sqrt(2)*$d;
  } else {
    print "BAD BOX OR BOX NOT DETECTED; USING DUMMY BOX SIZE IN GRO FILE";
    printf GRO "%11.5f%11.5f%11.5f\n", 1.0, 1.0, 1.0;
    #DLM made above change (and commented out below) 10/27/2007
    #die "Box type not supported."
  }

  close(GRO);
}

############    write .g96 file 

printf G96 "TITLE\n";
printf G96 "$outname.gro created by rdparm2gmx.pl $date\n $nato\n";
printf G96 "END\n";
printf G96 "POSITION\n";
for(my $i=1;$i<=$nato;$i++) {
  printf G96 "%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",$resnum{$i},$resname{$i},$atomname{$i},$i,$coordinates{"${i}x"},$coordinates{"${i}y"},$coordinates{"${i}z"};
}
printf G96 "END\n";
printf G96 "BOX\n";
# rectangular box
if($box_angles{1} == "90.0000000" && $box_angles{2} == "90.0000000" && $box_angles{3} == "90.0000000"){
  print "Rectangular box.";
  printf G96 "%15.9f%15.9f%15.9f\n", $box{x}, $box{y}, $box{z};
} elsif($box_angles{1} == "109.4712190" && $box_angles{2} == "109.4712190" && $box_angles{3} == "109.4712190") {
  print "Octahedral box.";
  # TODO: Check to make sure all box dimensions are equal.
  my $d = 2.0 * $box{x} + 0.0;
#  printf G96 "%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n", $d, (2./3.)*sqrt(2)*$d, (1./3.)*sqrt(6)*$d, 0.0, 0.0, (1.0/3.0)*$d, 0.0, -(1./3.)*$d, (1./3.)*sqrt(2)*$d;
  printf G96 "%15.9f%15.9f%15.9f\n", $d, (2./3.)*sqrt(2)*$d, (1./3.)*sqrt(6)*$d
#  my @ucell = ();
#  my @box = ($box{x}+0.0, $box{y}+0.0, $box{z}+0.0, $box_angles{1}+0.0, $box_angles{2}+0.0, $box_angles{3}+0.0);
#  my $DEGRAD = 0.0174532925;
#  $ucell[0] = $box[0]; 
#  $ucell[1] = 0.0; 
#  $ucell[2] = 0.0; 
#  $ucell[3] = $box[1]*cos($DEGRAD*$box[5]); 
#  $ucell[4] = $box[1]*sin($DEGRAD*$box[5]); 
#  $ucell[5] = 0.0; 
#  $ucell[6] = $box[2]*cos($DEGRAD*$box[4]); 
#  $ucell[7] = ($box[1]*$box[2]*cos($DEGRAD*$box[3]) - $ucell[6]*$ucell[3]) / $ucell[4]; 
#  $ucell[8] = sqrt($box[2]*$box[2] - $ucell[6]*$ucell[6] - $ucell[7]*$ucell[7]); 
#  printf G96 "%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n", @ucell;
} else {
  print "BAD BOX OR BOX NOT DETECTED; USING DUMMY BOX SIZE IN GRO FILE";
  printf G96 "%15.9f%15.9f%15.9f\n", 1.0, 1.0, 1.0;
  #DLM made above change (and commented out below) 10/27/2007
  #die "Box type not supported."
}
printf G96 "END\n";
close G96;

############	balance charges if not net-neutral by substracting from atom in protein with largest absolute charge ########
my $qtot = 0; # total charge
my $qmax = 0; # largest absolute charge
my $atmax = -1; # index of atom with largest absolute charge
for(my $i=1;$i<=$nato;$i++) {
  # accumulate total charge
  $qtot+=$charge{$i};

  if($resname{$i} ne "WAT" && $resname{$i} ne "Na+" && $resname{$i} ne "Cl-" && $resname{$i} ne "MG2") {
    # get absolute charge of this atom
    my $abch = abs($charge{$i});

    # keep track of solute atom with largest absolute charge
    if($abch > $qmax){ $qmax = $abch; $atmax = $i; } 
  }
}
print "\nQtot = $qtot\nQmax = $qmax at atom $atmax\n\n";

# If not integer charge, remove discrepancy from atom with largest net charge.
#Find nearest integer charge
my $intchg=int($qtot+0.5*($qtot <=> 0));

if(abs($qtot)-abs($intchg) != 0){
  #Difference in charges
  my $chgdiff=$qtot-$intchg;
  print "removing $chgdiff from charge on atom $atmax\n\n";
  $charge{$atmax}-=$chgdiff;
}

#==============================================================================
# Write parameter section of TOP file.
#==============================================================================

# make a list of which nonbonded atomtype each gaff atomtype is supposed to correspond to
print STDOUT "Constructing list of unique nonbonded atomtypes...\n";
my %lj_atomtypes_hash = ();
foreach my $atom_index ( keys %gaff ) {
  my $atom_type = $gaff{$atom_index};
  my $atom_type_index = $atom_type_indices{$atom_index};
#print "$atom_index $atom_type $atom_type_index\n";
  $lj_atomtypes_hash{$atom_type} = $atom_type_index;  
}
my @unique_atomtypes_list = keys %lj_atomtypes_hash;
printf STDOUT "%d unique atomtypes found.\n", $#unique_atomtypes_list + 1;

#my %all_atomtypes_hash = ();
#foreach my $atomtype ( values %gaff ) {
#  $all_atomtypes_hash{$atomtype} = 1;
#}
#
#my @atomtypes_with_parameters = keys %rad;
#my %all_atomtypes_hash = ();
#foreach my $atomtype ( values %gaff ) {
#  $all_atomtypes_hash{$atomtype} = 1;
#}
#my @all_atomtypes = keys %all_atomtypes_hash;
#
#print "all atomtypes: @all_atomtypes\n";
#print "atomtypes with parameters: @atomtypes_with_parameters\n";
#print "missing atomtypes: \n";
#foreach my $atomtype ( @all_atomtypes ) {
#    if(! defined $rad{$atomtype} ) {
#       print "$atomtype\n";
#    }
#}

# Write atom types.
print STDOUT "Writing atom types...\n";
print TOP "
[ atomtypes ]
;name  bond_type    mass    charge   ptype          sigma      epsilon\n";

#foreach my $atomtype ( keys %rad ) {
#  my $radg = $rad{$atomtype};
#  my $epsg = $eps{$atomtype};
#  # write zero default charges, as they will be overridden by individual atom entries
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", $atomtype, $atomtype, $radg, $epsg;
#}

foreach my $atom_type ( keys %lj_atomtypes_hash ) {
  my $atom_type_index = $lj_atomtypes_hash{$atom_type};
#print "$atom_type $atom_type_index\n";
  my $radg = $rad{$atom_type_index};
  my $epsg = $eps{$atom_type_index};
  # write zero default charges, as they will be overridden by individual atom entries
  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", $atom_type, $atom_type, $radg, $epsg;
}


# FIX: Add in missing atomtypes.
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "N3", "N3", $rad{"N"}, $eps{"N"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "O", "O", $rad{"O2"}, $eps{"O2"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "CB", "CB", $rad{"C"}, $eps{"C"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "N2", "N2", $rad{"N"}, $eps{"N"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "CA", "CA", $rad{"C"}, $eps{"C"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "NA", "NA", $rad{"N"}, $eps{"N"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "CN", "CN", $rad{"C"}, $eps{"C"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "C*", "C*", $rad{"C"}, $eps{"C"};
#  printf TOP "%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", "CW", "CW", $rad{"C"}, $eps{"C"};

# Write default bond types.
#Adding "if nwat" here, DLM, 1-10-06
if($nwat > 0) {
  print TOP "
[ bondtypes ]
; i    j      func       b0          kb
OW    HW         1    0.09572   462750.4 ; TIP3P water
HW    HW         1    0.15136   462750.4 ; TIP3P water
";

  # Write default angle types.
  print TOP "
[ angletypes ]
;  i    j    k  func       th0       cth
HW  OW  HW           1   104.520    836.800 ; TIP3P water
HW  HW  OW           1   127.740      0.000 ; (found in crystallographic water with 3 bonds)
";
}

  #==============================================================================
  # Write solute section of TOP file.
  #==============================================================================

print TOP "\n[ moleculetype ]
; Name            nrexcl
solute             3\n";

# Write solute atom records.
print TOP "\n[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB\n";
for(my $i=1; $i <= $nato; $i++) {
  if( is_solute($resname{$i}) ) {
    printf TOP "%6d %10s %6d %6s %6s %6d %10.5f %10f\n", $i, $gaff{$i}, $resnum{$i}, $resname{$i}, $atomname{$i}, $i, $charge{$i}, $mass{$i};
  }
}

# Write bond records.
print TOP "\n[ bonds ]
;  ai    aj funct  r  k\n";
for(my $bondindex = 1; $bondindex <= $nbon; $bondindex++) {
  my $i = $atom1{$bondindex};
  my $j = $atom2{$bondindex};
  my $Keq = $bstr{$bondindex};
  my $Req = $leng{$bondindex};

  if( is_solute($resname{$i}) && is_solute($resname{$j}) ) {
    printf TOP "%5d%6d%6d%12.4e%12.4e\n", $i, $j, 1, $Req, $Keq;
  }
}

# Write pair records.
# JDC: What do these do???
print TOP "\n[ pairs ]
;  ai    aj funct\n";
for(my $pairindex = 1; $pairindex <= $npai; $pairindex++) {
  my $i = $pair1{$pairindex};
  my $j = $pair2{$pairindex};

  if( is_solute($resname{$i}) && is_solute($resname{$j}) ) {
    printf TOP "%6d %6d %6d\n", $i, $j, 1;     
  }
}

print TOP "\n[ angles ]
;  ai    aj    ak funct  theta   cth\n";
for(my $angleindex = 1; $angleindex <= $nang; $angleindex++) {

  my $i = $atmo1{$angleindex};
  my $j = $atmo2{$angleindex};
  my $k = $atmo3{$angleindex};

  if( is_solute($resname{$i}) && is_solute($resname{$j}) && is_solute($resname{$k}) ) {
    printf TOP "%5d%6d%6d%6d%12.4e%12.4e\n",$i,$j,$k,1,$ang{$angleindex},$str{$angleindex};
  }
}

# using rdparm file mixes imps with props, so all torsions are converted to RB's #
# print TOP "\n[ dihedrals ]
# ;  ai    aj    ak    al funct  phi  k  mult\n";
# for($i=1;$i<=$nimp;$i++){
#   printf TOP "    %-4s%-4s%-4s%-4s%3d%12.4e%12.4e%12.4e\t;\n",$iat{$i,1},$iat{$i,2},$iat{$i,3},$iat{$i,4},1,$impphase{$i},$impkd{$i},$imppn{$i};
# }


# Write pair records.
#print TOP "\n[ exclusions ]
#;  ai    aj funct\n";
#for(my $pairindex = 1; $pairindex <= $npai; $pairindex++) {
#  my $i = $pair1{$pairindex};
#  my $j = $pair2{$pairindex};
#
#  if( is_solute($resname{$i}) && is_solute($resname{$j}) ) {
#    printf TOP "%6d %6d %6d\n", $i, $j, 1;     
#  }
#}

############ 	process proper dihedrals	################
# looks for multiple dihed's for same quartet #
sub same_line{
  my $same = 0;
  my $current = shift(@_);
  my $add = shift(@_);
  my $new = $current + $add;
  for(my $d=1; $d < 5; $d++) {
    if($at{$current,$d} eq $at{$new,$d}){ $same++; }
  }
  return($same);
}

my $d1type = 1;
my $func = 3;
for(my $i=1; $i<=$numpropers; $i++) { 
  # initialize the force constants for summing #
  my @V;
  my @C;
  for(my $j=0; $j<6; $j++) {
    $V[$j] = 0;
    $C[$j] = 0;
  }
  # print header if needed #
  if($d1type) { 
    print TOP "\n[ dihedrals ]\n\;i  j   k  l	 func	C0  ...  C5\n"; 
    $d1type = 0; 
  }

  # test lines that follow line $i #
  my $numijkl = 1;
  for(my $j=1; ($j<5) && (($i+$j) <= $numpropers); $j++) {
    my $test = same_line($i,$j);
    if($test == 4) { 
      $numijkl++; 
    }
  }
  
  # get all force constants for each line of a dihedral #
  my $lines = $i -1 +$numijkl;
  for(my $j=$i;$j<=$lines;$j++){
    my $period = abs($pn{$j});
    if($pk{$j}>0) { 
      $V[$period] = 2*$pk{$j}*$cal/$idivf{$j}; 
    }

    # assign V values to C values as predefined #
    if($period==1){
      $C[0]+=0.5*$V[$period];
      if($phase{$j}==0){
        $C[1]-=0.5*$V[$period];
      }else{
        $C[1]+=0.5*$V[$period];
      }
    }elsif($period==2){
      if(($phase{$j}==180)||($phase{$j}==3.14)){
        $C[0]+=$V[$period];
        $C[2]-=$V[$period];
      }else{
        $C[2]+=$V[$period];
      }
    }elsif($period==3){
      $C[0]+=0.5*$V[$period];
      if($phase{$j}==0){
        $C[1]+=1.5*$V[$period];
        $C[3]-=2*$V[$period];
      }else{
        $C[1]-=1.5*$V[$period];
        $C[3]+=2*$V[$period];
      }
    }elsif($period==4){
      if(($phase{$j}==180)||($phase{$j}==3.14)){
        $C[2]+=4*$V[$period];
        $C[4]-=4*$V[$period];
      }else{
        $C[0]+=$V[$period];
        $C[2]-=4*$V[$period];
        $C[4]+=4*$V[$period];
      }
    }
  }
  printf TOP "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\t;\n",$at{$i,1},$at{$i,2},$at{$i,3},$at{$i,4},$func,$C[0],$C[1],$C[2],$C[3],$C[4],$C[5];
  
  # skip the lines just analyzed #
  $i+=$numijkl-1;
}

#==============================================================================
# Write ion section of TOP file.
#==============================================================================

if($nsodium > 0) {
print TOP "
[ moleculetype ]
; molname       nrexcl
Na+             1

[ atoms ]
; id    at type res nr  residu name     at name  cg nr  charge   mass
1       IP           1          Na+         Na+       1      1   22.9898
"; }

if($nmagnesium > 0) {
print TOP "
[ moleculetype ]
; molname       nrexcl
MG2             1

[ atoms ]
; id    at type res nr  residu name     at name  cg nr  charge   mass
1       MG              1         MG2      MG     1      +2   24.30000
";}

if($nchloride > 0) {
print TOP "
[ moleculetype ]
; molname       nrexcl
Cl-             1

[ atoms ]
; id    at type res nr  residu name     at name  cg nr  charge   mass
1       IM              1         Cl-       Cl-      1      -1   35.45300
";}


#==============================================================================
# Write solvent section of TOP file.
#==============================================================================

#Adding if nwat here, DLM, 1-9-06
if($nwat > 0) {
  print TOP "

[ moleculetype ]
; molname       nrexcl
WAT             1

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.834   16.00000
     2     HW      1     WAT    H1      1      0.417    1.00800
     3     HW      1     WAT    H2      1      0.417    1.00800

[ settles ]
; OW    funct   doh     dhh
1       1       0.09572 0.15139

[ exclusions ]
1       2       3
2       1       3
3       1       2
\n";
}#end DLM

#==============================================================================
# Write system section of TOP file.
#==============================================================================

print TOP "

[ system ]
$nato system

[ molecules ]\n";
print TOP "; Compound        nmols\n";
print TOP "solute            1\n";
print TOP "Na+               $nsodium\n" if($nsodium > 0);
print TOP "MG2               $nmagnesium\n" if($nmagnesium > 0);
print TOP "Cl-               $nchloride\n" if($nchloride > 0);
print TOP "WAT               $nwat\n" if ($nwat > 0);

close(TOP);

# Clean up by removing intermediate files.
{
  my $command = "rm -f $rdparm_input_filename $rdparm_output_filename $pdb_filename";
  print "Removing temporary files...\n$command\n" if($debug_flag);
  system($command);
}

print "Done.\n";


