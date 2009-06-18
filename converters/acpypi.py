#!/usr/bin/env python
"""
    This code is released under GNU General Public License V3.

          <<<  NO WARRANTY AT ALL!!!  >>>

    It was inspired by:

    - amb2gmx.pl (Eric Sorin, David Mobley and John Chodera)
      and depends on Antechamber and Openbabel

    - YASARA Autosmiles:
      http://www.yasara.org/autosmiles.htm (Elmar Krieger)

    - topolbuild (Bruce Ray)

    - xplo2d (G.J. Kleywegt)

    For Antechamber, please cite:
    1.  Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and
        bond type perception in molecular mechanical calculations". Journal of
        Molecular Graphics and Modelling , 25, 2006, 247260.
    2.  Wang, J., Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A.
        "Development and testing of a general AMBER force field". Journal of
        Computational Chemistry, 25, 2004, 1157-1174.

    If you use this code, I am glad if you cite:

    BATISTA, P. R.; WILTER, A.; DURHAM, E. H. A. B.; PASCUTTI, P. G. Molecular
    Dynamics Simulations Applied to the Study of Subtypes of HIV-1 Protease.
    Cell Biochemistry and Biophysics, 44, 395-404, 2006.

    Alan Wilter S. da Silva, D.Sc. - CCPN Research Associate
    Department of Biochemistry, University of Cambridge.
    80 Tennis Court Road, Cambridge CB2 1GA, UK.
    >>http://www.bio.cam.ac.uk/~awd28<<

    alanwilter _at_ gmail _dot_ com
"""
svnId = '$Id: acpypi.py 234 2009-05-15 13:48:23Z alanwilter $'
try: svnRev,svnDate,svnTime = svnId.split()[2:5]
except: svnRev,svnDate,svnTime = 0,0,0
tag ="%s %s Rev: %s" % (svnDate, svnTime, svnRev)

lineHeader = \
'''
| ACPYPI: AnteChamber PYthon Parser Interface v. %s (c) 2008-2009 AWSdS |
''' % tag
frameLine = (len(lineHeader) -2) * '='
header = '%s%s%s' % (frameLine, lineHeader, frameLine)

#    TODO:
#        Howto Charmm and Amber with NAMD
#        Howto build topology for a modified amino acid
#        CYANA topology files

from commands import getoutput
from datetime import datetime
from shutil import copy2
from shutil import rmtree
import traceback
import signal
import time
import getopt
import math
import os
import cPickle as pickle
import string
import sys
import subprocess as sub

# List of Topology Formats created by acpypi so far:
outTopols = ['gmx', 'cns', 'charmm']

leapGaffFile = 'leaprc.gaff'
leapAmberFile = 'leaprc.ff99SB' #'leaprc.ff03.r1' #'oldff/leaprc.ff99'

cal = 4.184
Pi = 3.141593
qConv = 18.2223
radPi = 57.295780 # 180/Pi
maxDist = 3.0
minDist = 0.5
maxDist2 = maxDist**2 #squared Ang.
minDist2 = minDist**2 #squared Ang.

dictAmbAtomType2AmbGmxCode = \
{'BR':'1','C':'2','CA':'3','CB':'4','CC':'5','CK':'6','CM':'7','CN':'8','CQ':'9',
 'CR':'10','CT':'11','CV':'12','CW':'13','C*':'14','Ca':'15','F':'16','H':'17',
 'HC':'18','H1':'19','H2':'20','H3':'21','HA':'22','H4':'23','H5':'24','HO':'25',
 'HS':'26','HW':'27','HP':'28','I':'29','Cl':'30','Na':'31','IB':'32','Mg':'33',
 'N':'34','NA':'35','NB':'36','NC':'37','N2':'38','N3':'39','N*':'40','O':'41',
 'OW':'42','OH':'43','OS':'44','O2':'45','P':'46','S':'47','SH':'48','CU':'49',
 'FE':'50','K':'51','Rb':'52','Cs':'53','Li':'56','Zn':'57','Sr':'58','Ba':'59',
 'MCH3A':'MCH3A','MCH3B':'MCH3B','MNH2':'MNH2','MNH3':'MNH3','MW':'MW'}

dictOplsAtomType2OplsGmxCode = \
{'Ac3+':['697'],'Am3+':['699'],'Ar':['Ar','097'],'Ba2+':['414'],
 'Br':['722','730'],'Br-':['402'],
 'CT':['064','076','122','135','136','137','138','139','148','149','152','157','158','159','161','173','174','175','181','182','183','184','206','207','208','209','210','211','212','213','214','215','216','217','218','219','220','223','224','225','229','230','242','243','244','256','257','258','259','273','274','275','276','291','292','293','294','297','305','306','307','308','331','371','373','375','391','396','421','431','443','448','453','455','458','461','468','476','482','484','486','490','491','492','498','499','505','515','516','645','670','671','672','673','674','675','676','677','678','679','680','681','701','725','747','748','755','756','757','758','762','764','765','766','774','775','776','782','783','903','904','905','906','907','908','912','913','914','915','942','943','944','945','951','957','959','960','961','962','963','964'],
 'CA':['053','145','147','166','199','221','228','260','263','266','302','312','315','317','336','351','362','380','457','460','463','472','488','521','522','523','528','532','533','538','539','551','582','590','591','592','593','604','605','606','607','608','609','610','611','612','625','644','647','648','649','650','651','652','714','716','718','720','724','727','729','731','735','736','737','738','739','742','752','768','916','917','918'],
 'C3':['007','010','036','039','063','065','067','068','069','070','080','088','090','092','096','106','107','109','126','132','415','418','425','429'],
 'C':['001','017','026','058','095','131','231','234','235','247','252','267','320','322','334','366','378','470','471','772','952'],
 'C2':['005','009','015','016','019','022','027','028','031','034','037','056','057','061','071','081','089','091','093','110'],
 'CT_2':['223B','224B','225B','246','283','284','285','292B','293B','295','298','299','906B','912B'],
 'CM':['141','142','143','227','323','324','337','338','381','382','517','518','708'],
 'CW':['508','514','543','552','561','567','575','583','588','637'],
 'CB':['050','349','350','364','365','501','595','623','624'],
 'CH':['006','008','014','025','029','030','060','073'],
 'CZ':['261','423','754','925','927','928','929','931'],
 'CO':['189','191','193','195','197','198'],
 'C_2':['232','233','277','280','465'],
 'CR':['506','509','558','572','634'],
 'CQ':['347','531','621','642'],
 'CV':['507','560','574','636'],
 'CY':['711','712','713','733'],
 'CS':['544','568','589'],
 'CK':['353','627'],'CN':['502','594'],'CP':['043','048'],'CU':['550','581'],
 'CT_3':['245','296'],'C=':['150','178'],'CD':['011','075'],
 'C4':['066'],'C7':['077'],'C8':['074'],'C9':['072'],'CX':['510'],
 'C!':['145B'],'C*':['500'],'C+':['700'],'C_3':['271'],
 'CC':['045'],'CF':['044'],'CG':['049'],'CT_4':['160'],
 'Ca2+':['412'],
 'Cl':['123','151','226','264'],
 'Cl-':['401','709'],
 'Cs+':['410'],'Cu2+':['Cu2+'],'Eu3+':['705'],
 'F':['164','719','721','726','728','786','956','965'],
 'F-':['400'],'Fe2+':['Fe2+'],'Gd3+':['706'],
 'HA':['146','316','318','389','524','525','526','529','534','535','536','540','541','546','547','554','555','556','563','564','565','569','570','576','577','578','584','585','586','597','598','599','600','601','602','613','614','615','616','617','618','619','629','630','631','638','639','640','643','653','654','655','656','715','717','740','741','746'],
 'HC':['140','144','153','156','165','176','185','190','192','194','196','255','279','282','329','330','332','344','372','374','376','392','416','419','422','426','430','432','444','449','454','456','459','462','469','477','483','485','487','702','710','759','763','777','778','779','784','911','926','930','950','958'],
 'H':['004','013','041','047','128','240','241','250','254','314','325','327','339','342','343','357','358','360','367','369','383','385','387','388','428','479','481','504','513','545','553','562','596','632','744','745','909','910'],
 'H3':['021','052','055','104','105','289','290','301','304','310','941','955'],
 'HO':['024','079','155','163','168','170','172','188','270','435'],
 'HS':['033','086','087','204','205'],
 'HW':['112','114','117','119','796'],
 'H4':['345','390'],'H5':['355','359'],
 'He':['130'],'I':['732'],'I-':['403'],'K+':['408'],'Kr':['098'],
 'LP':['433','797'],'La3+':['703'],'Li+':['404','406'],
 'MCH3A':['MCH3A'],'MCH3B':['MCH3B'],'MNH2':['MNH2'],'MNH3':['MNH3'],
 'MW':['MW','115'],'Mg2+':['411'],
 'NA':['040','046','319','321','333','354','361','377','379','503','512','542','548','557','587','628'],
 'NC':['311','335','346','348','363','520','527','530','537','603','620','622','641','646'],
 'N':['003','012','094','237','238','239','249','251','265','478','480','787'],
 'N3':['020','101','102','103','286','287','288','309','427','940','953'],
 'N2':['051','054','300','303','313','341','356','368','386','743'],
 'NB':['042','352','511','549','559','573','580','626','635'],
 'N*':['319B','333B','354B','377B'],
 'NT':['127','900','901','902'],
 'NZ':['262','424','750','753'],
 'NO':['760','767'],'NY':['749','751'],
 'Na+':['405','407'],'Nd3+':['704'],'Ne':['129'],
 'OS':['062','108','179','180','186','395','442','447','452','467','473','566','571','579','773'],
 'O':['002','059','236','248','253','326','328','340','370','384','771','788'],
 'OH':['023','078','154','162','167','169','171','187','268','420','434'],
 'O2':['018','125','272','394','441','446','451','954'],
 'OW':['111','113','116','118','795'],
 'O_2':['278','281','466'],
 'OY':['475','494','497'],
 'OL':['120'],'ON':['761'],'OU':['437'],'O_3':['269'],
 'P':['393','440','445','450','785'],
 'P+':['781'],'Rb+':['409'],
 'S':['035','038','084','085','124','202','203','222','633'],
 'SH':['032','082','083','200','201','417','734'],
 'SI':['SI'],'SY':['474'],'SY2':['493'],'SZ':['496'],'Sr2+':['413'],
 'Th4+':['698'],'U':['436'],'Xe':['099'],'Yb3+':['707'],'Zn2+':['Zn2+']}

# reverse dictOplsAtomType2OplsGmxCode
oplsCode2AtomTypeDict = {}
for k,v in dictOplsAtomType2OplsGmxCode.items():
    for code in v:
        oplsCode2AtomTypeDict[code] = k
#        if oplsCode2AtomTypeDict.has_key(code):
#            oplsCode2AtomTypeDict[code].append(k)
#        else:
#            oplsCode2AtomTypeDict[code] = [k]

# Cross dictAmbAtomType2AmbGmxCode with dictOplsAtomType2OplsGmxCode & add H1,HP,H2
dictAtomTypeAmb2OplsGmxCode = {'H1':['140','1.00800'],'HP':['140','1.00800'],'H2':['140','1.00800']}
dictOplsMass = {'SY2': ['32.06000'], 'Zn2+': ['65.37000'], 'CQ': ['12.01100'], 'CP': ['12.01100'], 'Nd3+': ['144.24000'], 'Br-': ['79.90400'], 'Cu2+': ['63.54600'], 'Br': ['79.90400'], 'H': ['1.00800'], 'P': ['30.97376'], 'Sr2+': ['87.62000'], 'ON': ['15.99940'], 'OL': ['0.00000'], 'OH': ['15.99940'], 'OY': ['15.99940'], 'OW': ['15.99940'], 'OU': ['15.99940'], 'OS': ['15.99940'], 'Am3+': ['243.06000'], 'HS': ['1.00800'], 'HW': ['1.00800'], 'HO': ['1.00800'], 'HC': ['1.00800'], 'HA': ['1.00800'], 'O2': ['15.99940'], 'Ca2+': ['40.08000'], 'Th4+': ['232.04000'], 'He': ['4.00260'], 'C': ['12.01100'], 'Cs+': ['132.90540'], 'O': ['15.99940'], 'Gd3+': ['157.25000'], 'S': ['32.06000'], 'P+': ['30.97376'], 'La3+': ['138.91000'], 'H3': ['1.00800'], 'H4': ['1.00800'], 'MNH2': ['0.00000'], 'MW': ['0.00000'], 'NB': ['14.00670'], 'K+': ['39.09830'], 'Ne': ['20.17970'], 'Rb+': ['85.46780'], 'C+': ['12.01100'], 'C*': ['12.01100'], 'NO': ['14.00670'], 'CT_4': ['12.01100'], 'NA': ['14.00670'], 'C!': ['12.01100'], 'NC': ['14.00670'], 'NZ': ['14.00670'], 'CT_2': ['12.01100'], 'CT_3': ['12.01100'], 'NY': ['14.00670'], 'C9': ['14.02700'], 'C8': ['13.01900'], 'C=': ['12.01100'], 'Yb3+': ['173.04000'], 'C3': ['15.03500', '12.01100'], 'C2': ['14.02700'], 'C7': ['12.01100'], 'C4': ['16.04300'], 'CK': ['12.01100'], 'Cl-': ['35.45300'], 'N*': ['14.00670'], 'CH': ['13.01900'], 'CO': ['12.01100'], 'CN': ['12.01100'], 'CM': ['12.01100'], 'F': ['18.99840'], 'CC': ['12.01100'], 'CB': ['12.01100'], 'CA': ['12.01100'], 'CG': ['12.01100'], 'CF': ['12.01100'], 'N': ['14.00670'], 'CZ': ['12.01100'], 'CY': ['12.01100'], 'CX': ['12.01100'], 'Ac3+': ['227.03000'], 'CS': ['12.01100'], 'CR': ['12.01100'], 'N2': ['14.00670'], 'N3': ['14.00670'], 'CW': ['12.01100'], 'CV': ['12.01100'], 'CU': ['12.01100'], 'CT': ['12.01100'], 'SZ': ['32.06000'], 'SY': ['32.06000'], 'Cl': ['35.45300'], 'NT': ['14.00670'], 'O_2': ['15.99940'], 'Xe': ['131.29300'], 'SI': ['28.08000'], 'SH': ['32.06000'], 'Eu3+': ['151.96000'], 'F-': ['18.99840'], 'MNH3': ['0.00000'], 'H5': ['1.00800'], 'C_3': ['12.01100'], 'C_2': ['12.01100'], 'I-': ['126.90450'], 'LP': ['0.00000'], 'I': ['126.90450'], 'Na+': ['22.98977'], 'Li+': ['6.94100'], 'U': ['0.00000'], 'MCH3A': ['0.00000'], 'MCH3B': ['0.00000'], 'CD': ['13.01900', '12.01100'], 'O_3': ['15.99940'], 'Kr': ['83.79800'], 'Fe2+': ['55.84700'], 'Ar': ['39.94800'], 'Mg2+': ['24.30500'], 'Ba2+': ['137.33000']}
for ambKey in dictAmbAtomType2AmbGmxCode.keys():
    if dictOplsAtomType2OplsGmxCode.has_key(ambKey):
        dictAtomTypeAmb2OplsGmxCode[ambKey] = dictOplsAtomType2OplsGmxCode[ambKey]+list(dictOplsMass[ambKey])

# learnt from 22 residues test.
dictAtomTypeAmb2OplsGmxCode = {'HS': ['204', '1.008'], 'HP': ['140', '1.008'], 'HO': ['155', '168', '1.008'], 'HC': ['140', '1.008'], 'HA': ['146', '1.008'], 'O2': ['272', '15.9994'], 'C*': ['500', '12.011'], 'NA': ['503', '512', '14.0067'], 'NB': ['511', '14.0067'], 'CB': ['501', '12.011'], 'C': ['235', '271', '12.011'], 'CN': ['502', '12.011'], 'CM': ['302', '12.011'], 'CC': ['507', '508', '510', '12.011'], 'H': ['240', '241', '290', '301', '304', '310', '504', '513', '1.008'], 'CA': ['145', '166', '12.011'], 'O': ['236', '15.9994'], 'N': ['237', '238', '239', '14.0067'], 'S': ['202', '32.06'], 'CR': ['506', '509', '12.011'], 'N2': ['300', '303', '14.0067'], 'N3': ['287', '309', '14.0067'], 'CW': ['508', '510', '514', '12.011'], 'CV': ['507', '12.011'], 'CT': ['135', '136', '137', '149', '157', '158', '206', '209', '210', '223B', '224B', '245', '246', '274', '283', '284', '285', '292', '292B', '293B', '296', '307', '308', '505', '12.011'], 'OH': ['154', '167', '15.9994'], 'H1': ['140', '1.008'], 'H4': ['146', '1.008'], 'H5': ['146', '1.008'], 'SH': ['200', '32.06']}

# learnt from 22 residues test.
dictAtomTypeGaff2OplsGmxCode = {'cc': ['500', '506', '507', '508', '514', '12.011'], 'ca': ['145', '166', '501', '502', '12.011'], 'h1': ['140', '1.008'], 'h4': ['146', '1.008'], 'h5': ['146', '1.008'], 'cz': ['302', '12.011'], 'c2': ['509', '510', '12.011'], 'nh': ['300', '303', '14.0067'], 'ha': ['146', '1.008'], 'na': ['503', '512', '14.0067'], 'nc': ['511', '14.0067'], 'nd': ['511', '14.0067'], 'hx': ['140', '1.008'], 'hs': ['204', '1.008'], 'hn': ['240', '241', '290', '301', '304', '310', '504', '513', '1.008'], 'ho': ['155', '168', '1.008'], 'c3': ['135', '136', '137', '149', '157', '158', '206', '209', '210', '223B', '224B', '245', '246', '274', '283', '284', '285', '292', '292B', '293B', '296', '307', '308', '505', '12.011'], 'hc': ['140', '1.008'], 'cd': ['500', '506', '507', '508', '514', '12.011'], 'c': ['235', '271', '12.011'], 'oh': ['154', '167', '15.9994'], 'ss': ['202', '32.06'], 'o': ['236', '272', '15.9994'], 'n': ['237', '238', '239', '14.0067'], 'sh': ['200', '32.06'], 'n4': ['287', '309', '14.0067']}

#draft
atomTypeAmber2oplsDict = {'HS': ['HS'], 'HP': ['HC'], 'HO': ['HO'], 'HC': ['HC'],
                   'HA': ['HA'], 'O2': ['O2'], 'C*': ['C*'], 'NA': ['NA'],
                   'NB': ['NB'], 'CB': ['CB'], 'CN': ['CN'], 'CV': ['CV'],
                   'CM': ['CA'], 'CA': ['CA'], 'CR': ['CR'], 'OH': ['OH'],
                   'H1': ['HC'], 'H4': ['HA'], 'N2': ['N2'], 'N3': ['N3'],
                   'H5': ['HA'], 'SH': ['SH'], 'N': ['N'], 'S': ['S'], 'O': ['O'],
                   'C': ['C', 'C_3'], 'CW': ['CW', 'CX'], 'H': ['H', 'H3'],
                   'CC': ['CX', 'CW', 'CV'], 'CT': ['CT', 'CT_2', 'CT_3']}

# draft
a2oD = {'amber99_2': ['opls_235', 'opls_271'],
        'amber99_3': ['opls_302', 'opls_145'],
        'amber99_5': ['opls_507', 'opls_508', 'opls_510'],
        'amber99_11': ['opls_209','opls_158','opls_283','opls_223B','opls_293B',
                       'opls_284','opls_292B','opls_274','opls_136','opls_135',
                       'opls_292','opls_157','opls_206','opls_137','opls_505',
                       'opls_224B','opls_307','opls_308','opls_210','opls_149'],
        'amber99_13': ['opls_514'],
        'amber99_14': ['opls_500'],
        'amber99_17': ['opls_504','opls_241','opls_240','opls_290','opls_301',
                       'opls_310','opls_304','opls_513'],
        'amber99_18': ['opls_140'],
        'amber99_19': ['opls_140'],
        'amber99_22': ['opls_146'],
        'amber99_23': ['opls_146'],
        'amber99_25': ['opls_155'],
        'amber99_26': ['opls_204'],
        'amber99_28': ['opls_140'],
        'amber99_34': ['opls_238', 'opls_239', 'opls_237'],
        'amber99_35': ['opls_512', 'opls_503'],
        'amber99_36': ['opls_511'],
        'amber99_38': ['opls_300', 'opls_303'],
        'amber99_39': ['opls_309', 'opls_287'],
        'amber99_41': ['opls_236'],
        'amber99_43': ['opls_154'],
        'amber99_45': ['opls_272'],
        'amber99_47': ['opls_202'],
        'amber99_48': ['opls_200'],
 }

global pid
pid = 0

head = "%s created by acpypi on %s\n"

date = datetime.now().ctime()

USAGE = \
"""
    acpypi -i _file_ [-c _string_] [-n _int_] [-m _int_] [-a _string_] [-f] etc. or
    acpypi -p _prmtop_ -x _inpcrd_ [-d]
    -i    input file name with either extension '.pdb' or '.mol2' (mandatory if -p and -x not set)
    -p    amber prmtop file name (always used with -x)
    -x    amber inpcrd file name (always used with -p)
    -c    charge method: gas, bcc (default), user (user's charges in mol2 file)
    -n    net molecular charge (int), for gas default is 0
    -m    multiplicity (2S+1), default is 1
    -a    atom type, can be gaff or amber (AMBER99SB), default is gaff
    -f    force topologies recalculation anew
    -d    for debugging purposes, keep any temporary file created
    -o    output topologies: all (default), gmx, cns
    -t    write CNS topology with allhdg-like parameters (experimental)
    -e    engine: tleap (default) or sleap (not fully matured)
    -b    a basename for the project (folder and output files)
    -s    max time (in sec) tolerance for mopac, default is 10 hours

    output: assuming 'root' is the basename of either the top input file,
            the 3-letter residue name or user defined (-b option)
    root_bcc_gaff.mol2:  final mol2 file with 'bcc' charges and 'gaff' atom type
    root_AC.inpcrd    :  coord file for AMBER
    root_AC.prmtop    :  topology and parameter file for AMBER
    root_AC.lib       :  residue library file for AMBER
    root_AC.frcmod    :  modified force field parameters
    root_GMX.gro      :  coord file for GROMACS
    root_GMX.top      :  topology file for GROMACS
    root_GMX.itp      :  molecule unit topology and parameter file for GROMACS
    root_GMX_OPLS.itp :  OPLS/AA mol unit topol & par file for GROMACS (experimental!)
    em.mdp, md.mdp    :  run parameters file for GROMACS
    root_NEW.pdb      :  final pdb file generated by ACPYPI
    root_CNS.top      :  topology file for CNS/XPLOR
    root_CNS.par      :  parameter file for CNS/XPLOR
    root_CNS.inp      :  run parameters file for CNS/XPLOR
    root_CHARMM.rtf   :  topology file for CHARMM
    root_CHARMM.prm   :  parameter file for CHARMM
    root_CHARMM.inp   :  run parameters file for CHARMM
"""

SLEAP_TEMPLATE = \
"""
source %(leapAmberFile)s
source %(leapGaffFile)s
set default fastbld on
#set default disulfide auto
%(res)s = loadpdb %(base)s.pdb
#check %(res)s
saveamberparm %(res)s %(acBase)s.prmtop %(acBase)s.inpcrd
saveoff %(res)s %(acBase)s.lib
quit
"""

TLEAP_TEMPLATE = \
"""
verbosity 1
source %(leapAmberFile)s
source %(leapGaffFile)s
mods = loadamberparams %(acBase)s.frcmod
%(res)s = loadmol2 %(acMol2FileName)s
check %(res)s
saveamberparm %(res)s %(acBase)s.prmtop %(acBase)s.inpcrd
saveoff %(res)s %(acBase)s.lib
quit
"""

def invalidArgs(text = None, usage = True):

    if usage: print USAGE
    if text:
        print 'ERROR: ' + text
    sys.exit(1)

verNum = string.split(sys.version)[0]
version = string.split(verNum, ".")
if map(int, version) < [2, 4, 0]:
    invalidArgs(text = "Python version %s\n       Sorry, you need python 2.4 or higher" % verNum, usage = False)

try:
    set()
except NameError:
    from sets import Set as set

def parseArgs(args):

    for _i in range(args.count('')):
        args.remove('')

    amb2gmx = False

    options = 'hftdi:c:n:m:o:a:e:x:p:b:s:'

    ctList = ['gas', 'bcc', 'user']
    atList = ['gaff', 'amber'] #, 'bcc', 'sybyl']
    tpList = ['all'] + outTopols
    enList = ['sleap', 'tleap']

    if '-h' in args:
        invalidArgs()
    try:
        opt_list, args = getopt.getopt(args, options) #, long_options)
    except Exception, msg:
        invalidArgs("FAILED: %s" % msg)

    d = {}

    for key, value in opt_list:
        if key in d:
            invalidArgs("duplicated parameters: %s" % key)

        if value == '':
            value = None

        d[key] = value

    if not d and not args:
        invalidArgs("missing parameters")

    if '-c' in d.keys():
        if d['-c'] not in ctList:
            invalidArgs("'%s' not valid, has to be one of %s" % (d['-c'], str(ctList)))

    if '-a' in d.keys():
        if d['-a'] not in atList:
            invalidArgs("'%s' not valid, has to be one of %s" % (d['-a'], str(atList)))

    if '-o' in d.keys():
        if d['-o'] not in tpList:
            invalidArgs("'%s' not valid, has to be one of %s" % (d['-o'], str(tpList)))

    if '-e' in d.keys():
        if d['-e'] not in enList:
            invalidArgs("'%s' not valid, has to be one of %s" % (d['-e'], str(enList)))

    if not '-i' in d.keys():
        amb2gmx = True
        if not '-p' in d.keys() or not '-x' in d.keys():
            invalidArgs("missing input files")

    if '-i' in d.keys() and '-p' in d.keys():
        invalidArgs("either '-i' or ('-p', '-x'), but not both")
    elif '-i' in d.keys() and '-x' in d.keys():
        invalidArgs("either '-i' or ('-p', '-x'), but not both")

    if amb2gmx:
        if not os.path.exists(d['-p']) or not os.path.exists(d['-x']):
            invalidArgs("input file prmtop and/or inpcrd don't exist")
    elif not os.path.exists(d['-i']):
        invalidArgs("input file '%s' doesn't exist" % d['-i'])

    if args:
        invalidArgs("argument not recognised: %s" %str(args))

    d['amb2gmx'] = amb2gmx

    return d

def elapsedTime(seconds, suffixes=['y','w','d','h','m','s'], add_s=False, separator=' '):
    """
    Takes an amount of seconds and turns it into a human-readable amount of time.
    """
    # the formatted time string to be returned
    time = []

    # the pieces of time to iterate over (days, hours, minutes, etc)
    # - the first piece in each tuple is the suffix (d, h, w)
    # - the second piece is the length in seconds (a day is 60s * 60m * 24h)
    parts = [(suffixes[0], 60 * 60 * 24 * 7 * 52),
          (suffixes[1], 60 * 60 * 24 * 7),
          (suffixes[2], 60 * 60 * 24),
          (suffixes[3], 60 * 60),
          (suffixes[4], 60),
          (suffixes[5], 1)]

    # for each time piece, grab the value and remaining seconds, and add it to
    # the time string
    for suffix, length in parts:
        value = seconds / length
        if value > 0:
            seconds = seconds % length
            time.append('%s%s' % (str(value),
                           (suffix, (suffix, suffix + 's')[value > 1])[add_s]))
        if seconds < 1:
            break

    return separator.join(time)

def splitBlock(dat):
    '''split a amber parm dat file in blocks
       0 = mass, 1 = extra + bond, 2 = angle, 3 = dihedral, 4 = improp, 5 = hbond
       6 = equiv nbon, 7 = nbon, 8 = END, 9 = etc.
    '''
    dict = {}
    count = 0
    for line in dat:
        line = line.rstrip()
        if dict.has_key(count):
            dict[count].append(line)
        else:
            dict[count] = [line]
        if not line:
            count += 1
    return dict

def getParCode(line):
    key = line.replace(' -','-').replace('- ','-').split()[0]
    return key

def parseFrcmod(lista):
    heads = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'HBON', 'NONB']
    dict = {}
    for line in lista[1:]:
        line = line.strip()
        if line[:4] in heads:
            head = line[:4]
            dict[head] = []
            dd = {}
            continue
        elif line:
            key = line.replace(' -','-').replace('- ','-').split()[0]
            if dd.has_key(key):
                if not dd[key].count(line):
                    dd[key].append(line)
            else:
                dd[key] = [line]
            dict[head]= dd
    for k in dict.keys():
        if not dict[k]: dict.pop(k)
    return dict

def parmMerge(fdat1, fdat2, frcmod = False):
    '''merge two amber parm dat/frcmod files and save in /tmp'''
    name1 = os.path.basename(fdat1).split('.dat')[0]
    if frcmod:
        name2 = os.path.basename(fdat2).split('.')[1]
    else:
        name2 = os.path.basename(fdat2).split('.dat')[0]
    mname = '/tmp/'+name1+name2+'.dat'
    mdatFile = open(mname,'w')
    mdat = ['merged %s %s' % (name1,name2)]

    #if os.path.exists(mname): return mname
    dat1 = splitBlock(file(fdat1).readlines())

    if frcmod:
        dHeads = {'MASS':0, 'BOND':1, 'ANGL':2, 'DIHE':3, 'IMPR':4, 'HBON':5, 'NONB':7}
        dat2 = parseFrcmod(file(fdat2).readlines()) # dict
        for k in dat2.keys():
            for parEntry in dat2[k].keys():
                idFirst = None
                for line in dat1[dHeads[k]][:]:
                    if line:
                        key = line.replace(' -','-').replace('- ','-').split()[0]
                        if key == parEntry:
                            if not idFirst: idFirst = dat1[dHeads[k]].index(line)
                            dat1[dHeads[k]].remove(line)
                rev = dat2[k][parEntry][:]
                rev.reverse()
                if idFirst == None: idFirst = 0
                for ll in rev:
                    dat1[dHeads[k]].insert(idFirst,ll)
        dat1[0][0] = mdat[0]
        for k in dat1.keys():
            for line in dat1[k]:
                mdatFile.write(line+'\n')
        return mname

    dat2 = splitBlock(file(fdat2).readlines())
    for k in dat1.keys()[:8]:
        if k == 0:
            lines = dat1[k][1:-1]+dat2[k][1:-1]+['']
            for line in lines:
                mdat.append(line)
        if k == 1:
            for i in dat1[k]:
                if '-' in i:
                    id1 = dat1[k].index(i)
                    break
            for j in dat2[k]:
                if '-' in j:
                    id2 = dat2[k].index(j)
                    break
            l1 = dat1[k][:id1]
            l2 = dat2[k][:id2]
            line = ''
            for item in l1 + l2:
                line += item.strip() + ' '
            mdat.append(line)
            lines = dat1[k][id1:-1]+dat2[k][id2:-1]+['']
            for line in lines:
                mdat.append(line)
        if k in [2,3,4,5,6]: # angles, p dih, imp dih
            lines = dat1[k][:-1]+dat2[k][:-1]+['']
            for line in lines:
                mdat.append(line)
        if k == 7:
            lines = dat1[k][:-1]+dat2[k][1:-1]+['']
            for line in lines:
                mdat.append(line)
    for k in dat1.keys()[8:]:
        for line in dat1[k]:
            mdat.append(line)
    for k in dat2.keys()[9:]:
        for line in dat2[k]:
            mdat.append(line)
    for line in mdat:
        mdatFile.write(line+'\n')
    mdatFile.close()

    return mname

class AbstractTopol:
    """
        Super class to build topologies
    """
    def __init__(self):
        if self.__class__ is AbstractTopol:
            raise TypeError("Attempt to create istance of abstract class AbstractTopol")

    def printDebug(self, text = ''):
        if self.debug:
            print 'DEBUG: %s' % text

    def printWarn(self, text = ''):
        print 'WARNING: %s' % text

    def printError(self, text = ''):
        print 'ERROR: %s' % text

    def printMess(self, text = ''):
        print '==> %s' % text

    def printQuoted(self, text = ''):
        print 10 * '+' + 'start_quote' + 59 * '+'
        print text
        print 10 * '+' + 'end_quote' + 61 * '+'

    def guessCharge(self):
        """
            Guess the charge of a system based on antechamber
            Returns None in case of error
        """
        done = False
        error = False
        charge = self.chargeVal
        localDir = os.path.abspath('.')
        if not os.path.exists(self.tmpDir):
            os.mkdir(self.tmpDir)
        if not os.path.exists(os.path.join(self.tmpDir, self.inputFile)):
            copy2(self.absInputFile, self.tmpDir)
        os.chdir(self.tmpDir)

        if self.chargeType == 'user':
            if self.ext == '.mol2':
                self.printMess("Reading user's charges from mol2 file...")
                charge = self.readMol2TotalCharge(self.inputFile)
                done = True
            else:
                self.printWarn("cannot read charges from a PDB file")
                self.printWarn("using now 'bcc' method for charge")

        if not self.chargeVal and not done:
            self.printWarn("no charge value given, trying to guess one...")
            mol2FileForGuessCharge = self.inputFile
            if self.ext == ".pdb":
                cmd = '%s -ipdb %s -omol2 %s.mol2' % (self.babelExe, self.inputFile,
                                                      self.baseName)
                self.printDebug("guessCharge: " + cmd)
                out = getoutput(cmd)
                self.printDebug(out)
                mol2FileForGuessCharge = os.path.abspath(self.baseName+".mol2")

            cmd = '%s -i %s -fi mol2 -o tmp -fo mol2 -c gas -pf y' % \
                                                        (self.acExe, mol2FileForGuessCharge)

            if self.debug:
                self.printMess("Debugging...")
                cmd = cmd.replace('-pf y', '-pf n')
            print cmd

            log = getoutput(cmd)
            #print log
            m = log.splitlines()
            if len(m) > 5:
                try:
                    n = m[-1].split()
                    charge = float(n[14].replace('(','').replace(')',''))
                except:
                    error = True
            elif len(m) == 0:
                self.printMess("An old version of Antechamber? Still trying to get charge...")
                charge = self.readMol2TotalCharge('tmp')
            else:
                error = True

            if error:
                self.printError("guessCharge failed")
                os.chdir(localDir)
                self.printQuoted(log)
                self.printMess("Trying with net charge = 0 ...")
#                self.chargeVal = 0
                return None

        charge = int(charge)
        self.chargeVal = str(charge)
        self.printMess("... charge set to %i" % charge)
        os.chdir(localDir)

    def setResNameCheckCoords(self):
        """Set a 3 letter residue name
           and check coords duplication
        """
        exit = False
        localDir = os.path.abspath('.')
        if not os.path.exists(self.tmpDir):
            os.mkdir(self.tmpDir)
        #if not os.path.exists(os.path.join(tmpDir, self.inputFile)):
        copy2(self.absInputFile, self.tmpDir)
        os.chdir(self.tmpDir)

        if self.ext == '.pdb':
            tmpFile = open(self.inputFile, 'r')
        else:
            cmd = '%s -i %s -fi %s -o tmp -fo ac -pf y' % \
                                    (self.acExe, self.inputFile, self.ext[1:])
            self.printDebug(cmd)
            out = getoutput(cmd)
            if not out.isspace():
                self.printDebug(out)
            try:
                tmpFile = open('tmp', 'r')
            except:
                rmtree(self.tmpDir)
                raise

        tmpData = tmpFile.readlines()
        residues = set()
        coords = {}
        for line in tmpData:
            if 'ATOM  ' in line or 'HETATM' in line:
                residues.add(line[17:20])
                at = line[0:17]
                cs = line[30:54]
                if coords.has_key(cs):
                    coords[cs].append(at)
                else:
                    coords[cs] = [at]
        #self.printDebug(coords)

        if len(residues) > 1:
            self.printError("more than one residue detected '%s'" % str(residues))
            self.printError("verify your input file '%s'. Aborting ..." % self.inputFile)
            sys.exit(1)

        dups = ""
        short = ""
        long = ""
        longSet = set()
        id = 0
        items = coords.items()
        l = len(items)
        for item in items:
            id += 1
            if len(item[1]) > 1: # if True means atoms with same coordinates
                for i in item[1]:
                    dups += "%s %s\n" % (i, item[0])

#        for i in xrange(0,len(data),f):
#            fdata += (data[i:i+f])+' '

            for id2 in xrange(id,l):
                item2 = items[id2]
                c1 = map(float,[item[0][i:i+8] for i in xrange(0,24,8)])
                c2 = map(float,[item2[0][i:i+8] for i in xrange(0,24,8)])
                dist2 = self.distance(c1,c2)
                if dist2 < minDist2:
                    dist = math.sqrt(dist2)
                    short += "%8.5f       %s %s\n" % (dist, item[1], item2[1])
                if dist2 < maxDist2: # and not longOK:
                    longSet.add(str(item[1]))
                    longSet.add(str(item2[1]))
            if str(item[1]) not in longSet:
                long += "%s\n" % item[1]

        if dups:
            self.printError("Atoms with same coordinates in '%s'!" % self.inputFile)
            self.printQuoted(dups[:-1])
            exit = True

        if short:
            self.printError("Atoms TOO close (< %s Ang.)" % minDist)
            self.printQuoted("Dist (Ang.)    Atoms\n" + short[:-1])
            exit = True

        if long:
            self.printError("Atoms TOO alone (> %s Ang.)" % maxDist)
            self.printQuoted(long[:-1])
            exit = True

        if exit:
            if self.force:
                self.printWarn("You chose to proceed anyway with '-f' option. GOOD LUCK!")
            else:
                self.printError("Use '-f' option if you want to proceed anyway. Aborting ...")
                rmtree(self.tmpDir)
                sys.exit(1)

        resname = list(residues)[0]
        newresname = resname

        if resname.isdigit() or 'E' in resname[1:3].upper() or 'ADD' in resname.upper():
            newresname = 'R' + resname
        if not resname.isalnum():
            newresname = 'MOL'
        if newresname != resname:
            self.printWarn("In %s.lib, residue name will be '%s' instead of '%s' elsewhere"
                           % (self.acBaseName, newresname, resname))

        self.resName = newresname

        os.chdir(localDir)
        self.printDebug("setResNameCheckCoords done")

    def distance(self, c1, c2):
        #print c1, c2
        dist2 = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[0] - c2[0])**2 + \
                (c1[2] - c2[2])**2
        #dist2 = math.sqrt(dist2)
        return dist2

    def readMol2TotalCharge(self, mol2File):
        """Reads the charges in given mol2 file and returns the total
        """
        charge = 0.0
        ll = []
        cmd = '%s -i %s -fi mol2 -o tmp -fo mol2 -c wc -cf tmp.crg -pf y' % \
                                                        (self.acExe, mol2File)
        if self.debug:
            self.printMess("Debugging...")
            cmd = cmd.replace('-pf y', '-pf n')

        self.printDebug(cmd)

        log = getoutput(cmd)

        if log.isspace():
            tmpFile = open('tmp.crg', 'r')
            tmpData = tmpFile.readlines()
            for line in tmpData:
                ll += line.split()
            charge = sum(map(float,ll))
        elif self.debug:
            self.printQuoted(log)

        self.printDebug("readMol2TotalCharge: " + str(charge))

        return charge

    def execAntechamber(self, chargeType = None, atomType = None):
        """
            To call Antechamber and execute it

Usage: antechamber -i   input file name
                   -fi  input file format
                   -o   output file name
                   -fo  output file format
                   -c   charge method
                   -cf  charge file name
                   -nc  net molecular charge (int)
                   -a   additional file name
                   -fa  additional file format
                   -ao  additional file operation
                        crd : only read in coordinate
                        crg: only read in charge
                        name  : only read in atom name
                        type  : only read in atom type
                        bond  : only read in bond type
                   -m   multiplicity (2S+1), default is 1
                   -rn  residue name, if not available in the input file, default is MOL
                   -rf  residue toplogy file name in prep input file, default is molecule.res
                   -ch  check file name in gaussian input file, default is molecule
                   -ek  empirical calculation (mopac or divcon) keyword in a pair of quotation marks
                   -gk  gaussian keyword in a pair of quotation marks
                   -df  use divcon flag, 1 - use divcon; 0 - use mopac (default)
                   -at  atom type, can be gaff, amber, bcc and sybyl, default is gaff
                   -du  check atom name duplications, can be yes(y) or no(n), default is yes
                   -j   atom type and bond type prediction index, default is 4
                        0    : no assignment
                        1    : atom type
                        2    : full  bond types
                        3    : part  bond types
                        4    : atom and full bond type
                        5    : atom and part bond type
                   -s   status information, can be 0 (brief), 1 (the default) and 2 (verbose)
                   -pf  remove the intermediate files: can be yes (y) and no (n), default is no
                   -i -o -fi and -fo must appear in command lines and the others are optional

                             List of the File Formats

                file format type  abbre. index | file format type abbre. index
                ---------------------------------------------------------------
                Antechamber        ac       1  | Sybyl Mol2         mol2    2
                PDB                pdb      3  | Modified PDB       mpdb    4
                AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6
                Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8
                Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10
                Gaussian Output    gout    11  | Mopac Output       mopout 12
                Alchemy            alc     13  | CSD                csd    14
                MDL                mdl     15  | Hyper              hin    16
                AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18
                Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20
                Divcon Input       divcrt  21  | Divcon Output      divout 22
                Charmm             charmm  23
                --------------------------------------------------------------

                AMBER restart file can only be read in as additional file.

                             List of the Charge Methods

                charge method     abbre.  index | charge method      abbre. index
                ----------------------------------------------------------------
                RESP               resp     1  |  AM1-BCC            bcc     2
                CM1                cm1      3  |  CM2                cm2     4
                ESP (Kollman)      esp      5  |  Mulliken           mul     6
                Gasteiger          gas      7  |  Read in charge     rc      8
                Write out charge   wc       9  |  Delete Charge      dc     10
                ----------------------------------------------------------------
        """
        global pid

        self.printMess("Executing Antechamber...")

        self.makeDir()

        ct = chargeType or self.chargeType
        at = atomType or self.atomType

        cmd = '%s -i %s -fi %s -o %s -fo mol2 -c %s -nc %s -m %s -s 2 -df 0 -at\
 %s -pf y' % (self.acExe, self.inputFile, self.ext[1:], self.acMol2FileName,
                     ct, self.chargeVal, self.multiplicity, at)

        if self.debug:
            self.printMess("Debugging...")
            cmd = cmd.replace('-pf y', '-pf n')

        self.printDebug(cmd)

        if os.path.exists(self.acMol2FileName) and not self.force:
            self.printMess("AC output file present... doing nothing")
        else:
            try: os.remove(self.acMol2FileName)
            except: pass
            p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE)
            pid = p.pid
            signal.signal(signal.SIGALRM, self.signal_handler)
            signal.alarm(self.timeTol)

            out = p.stdout.read()
            self.acLog = out

        if os.path.exists(self.acMol2FileName):
            self.printMess("* Antechamber OK *")
        else:
            self.printQuoted(self.acLog)
            return True

    def signal_handler(self, signum, frame):#, pid = 0):
        global pid
        pids = self.job_pids_family(pid)
        print pid, pids
        self.printMess("Timed out! Process %s killed, max exec time (%ss) exceeded" \
                        % (pids, self.timeTol))
        os.system('kill -15 %s' % pids)
        raise Exception("MOPAC taking too long to finish... aborting!")

    def job_pids_family(self,jpid):
        '''INTERNAL: Return all job processes (PIDs)'''
        pid = repr(jpid)
        dict_pids = {}
        pids = [pid]
        out = getoutput("ps -A -o uid,pid,ppid|grep %i" % os.getuid()).split('\n')
        for item in out:
            vec = item.split()
            dict_pids[vec[2]] = vec[1]
        while 1:
            try:
                pid = dict_pids[pid]
                pids.append(pid)
            except KeyError:
                break
        return string.join(pids)

    def delOutputFiles(self):
        delFiles = ['mopac.in', 'mopac.pdb', 'mopac.out', 'tleap.in','sleap.in',
                    'divcon.pdb', 'fixbo.log',
                    'addhs.log', 'ac_tmp_ot.mol2', 'frcmod.ac_tmp',
                    'fragment.mol2', self.tmpDir] #leap.log
        self.printMess("Removing temporary files...")
        for file in delFiles:
            file = os.path.join(self.absHomeDir, file)
            if os.path.exists(file):
                if os.path.isdir(file):
                    rmtree(file)
                else:
                    os.remove(file)

    def checkXyzAndTopFiles(self):
        fileXyz = self.acXyzFileName
        fileTop = self.acTopFileName
        if os.path.exists(fileXyz) and os.path.exists(fileTop):
            #self.acXyz = fileXyz
            #self.acTop = fileTop
            return True
        return False

    def execSleap(self):
        global pid
        self.makeDir()

        if self.ext == '.mol2':
            self.printWarn("Sleap doesn't work with mol2 files yet...")
            return True

        if self.chargeType != 'bcc':
            self.printWarn("Sleap works only with bcc charge method")
            return True

        if self.atomType != 'gaff':
            self.printWarn("Sleap works only with gaff atom type")
            return True

        sleapScpt = SLEAP_TEMPLATE % self.acParDict

        fp = open('sleap.in','w')
        fp.write(sleapScpt)
        fp.close()

        cmd = '%s -f sleap.in' % self.sleapExe

        if self.checkXyzAndTopFiles() and not self.force:
            self.printMess("Topologies files already present... doing nothing")
        else:
            try: os.remove(self.acTopFileName) ; os.remove(self.acXyzFileName)
            except: pass
            self.printMess("Executing Sleap...")
            self.printDebug(cmd)

            p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE)
            pid = p.pid
            signal.signal(signal.SIGALRM, self.signal_handler)
            signal.alarm(self.timeTol)

            out = p.stdout.read()

            self.sleapLog = out
            self.checkLeapLog(self.sleapLog)

            if self.checkXyzAndTopFiles():
                self.printMess(" * Sleap OK *")
            else:
                self.printQuoted(self.sleapLog)
                return True

    def execTleap(self):

        fail = False

        self.makeDir()

        if self.ext == ".pdb":
            self.printMess('... converting pdb input file to mol2 input file')
            if self.convertPdbToMol2():
                self.printError("convertPdbToMol2 failed")

        #print self.chargeVal

        if self.execAntechamber():
            self.printError("Antechamber failed")
            fail = True
            #sys.exit(1)

        if self.execParmchk():
            self.printError("Parmchk failed")
            fail = True
            #sys.exit(1)

        if fail:
            return True

        tleapScpt = TLEAP_TEMPLATE % self.acParDict

        fp = open('tleap.in','w')
        fp.write(tleapScpt)
        fp.close()

        cmd = '%s -f tleap.in' % self.tleapExe

        if self.checkXyzAndTopFiles() and not self.force:
            self.printMess("Topologies files already present... doing nothing")
        else:
            try: os.remove(self.acTopFileName) ; os.remove(self.acXyzFileName)
            except: pass
            self.printMess("Executing Tleap...")
            self.printDebug(cmd)
            self.tleapLog = getoutput(cmd)
            self.checkLeapLog(self.tleapLog)

        if self.checkXyzAndTopFiles():
            self.printMess("* Tleap OK *")
        else:
            self.printQuoted(self.tleapLog)
            return True

    def checkLeapLog(self, log):
        log = log.splitlines(True)
        check = ''
        block = False
        for line in log:
            #print "*"+line+"*"
            if "Checking '" in line:
                #check += line
                block = True
            if "Checking Unit." in line:
                block = False
            if block:
                check += line
        self.printQuoted(check[:-1])

    def locateDat(self, aFile):
        '''locate a file pertinent to $AMBERHOME/dat/leap/parm/'''
        amberhome = os.environ.get('AMBERHOME')
        if amberhome:
            aFileF = os.path.join(amberhome,'dat/leap/parm',aFile)
            if os.path.exists(aFileF): return aFileF
        aFileF = os.path.join(os.path.dirname(self.acExe),'../dat/leap/parm',aFile)
        if os.path.exists(aFileF): return aFileF
        return None

    def execParmchk(self):

        self.makeDir()
        cmd = '%s -i %s -f mol2 -o %s' % (self.parmchkExe, self.acMol2FileName,
                                          self.acFrcmodFileName)

        if self.atomType == 'amber':
            parm99file = self.locateDat('parm99.dat')
            gaffFile = self.locateDat('gaff.dat')
            frcmodff99SB = self.locateDat('frcmod.ff99SB')
            parm99gaffFile = parmMerge(parm99file,gaffFile)
            parm99SBgaffFile = parmMerge(parm99gaffFile,frcmodff99SB, frcmod = True)
            #parm99gaffFile = '/Users/alan/workspace/acpypi/ffamber_additions/parm99SBgaff.dat'
            cmd += ' -p %s' % parm99SBgaffFile

        self.parmchkLog = getoutput(cmd)

        self.printDebug(cmd)

        if os.path.exists(self.acFrcmodFileName):
            check = self.checkFrcmod()
            if check:
                self.printWarn("Couldn't determine all parameters:")
                self.printMess("From file '%s'\n" % self.acFrcmodFileName+check)
            else:
                self.printMess("* Parmchk OK *")
        else:
            self.printQuoted(self.parmchkLog)
            return True

    def checkFrcmod(self):
        check = ""
        frcmodContent = open(self.acFrcmodFileName, 'r').readlines()
        for line in frcmodContent:
            if "ATTN, need revision" in line:
                check += line
        return check

    def convertPdbToMol2(self):
        if self.ext == '.pdb':
            if self.execBabel():
                self.printError("convert pdb to mol2 via babel failed")
                return True

    def execBabel(self):

        self.makeDir()

        cmd = '%s -ipdb %s -omol2 %s.mol2' % (self.babelExe, self.inputFile,
                                                  self.baseName)
        self.printDebug(cmd)
        self.babelLog = getoutput(cmd)
        self.ext = '.mol2'
        self.inputFile = self.baseName+self.ext
        self.acParDict['ext'] = 'mol2'
        if os.path.exists(self.inputFile):
            self.printMess("* Babel OK *")
        else:
            self.printQuoted(self.babelLog)
            return True

    def makeDir(self):

        os.chdir(self.rootDir)
        self.absHomeDir = os.path.abspath(self.homeDir)
        if not os.path.exists(self.homeDir):
            os.mkdir(self.homeDir)
        os.chdir(self.homeDir)
        copy2(self.absInputFile, '.')

        return True

    def createACTopol(self):
        """
            If successful, Amber Top and Xyz files will be generated
        """
        #sleap = False
        if self.engine == 'sleap':
            if self.execSleap():
                self.printError("Sleap failed")
                self.printMess("... trying Tleap")
                if self.execTleap():
                    self.printError("Tleap failed")
        if self.engine == 'tleap':
            if self.execTleap():
                self.printError("Tleap failed")
                if self.extOld == '.pdb':
                    self.printMess("... trying Sleap")
                    self.ext = self.extOld
                    self.inputFile = self.baseName+self.ext
                    if self.execSleap():
                        self.printError("Sleap failed")
        if not self.debug:
            self.delOutputFiles()

    def createMolTopol(self):
        """
            Create molTop obj
        """
        self.topFileData = open(self.acTopFileName, 'r').readlines()
        self.molTopol = MolTopol(self)
        if self.outTopols:
            if 'cns' in self.outTopols:
                self.molTopol.writeCnsTopolFiles()
            if 'gmx' in self.outTopols:
                self.molTopol.writeGromacsTopolFiles()
            if 'charmm' in self.outTopols:
                self.writeCharmmTopolFiles()
        self.pickleSave()

    def pickleSave(self):
        """
            To restore:
                from acpypi import *
                import cPickle as pickle
                o = pickle.load(open('DDD.pkl'))
                NB: It fails to restore with ipython in Mac (Linux OK)
        """
        pklFile = self.baseName+".pkl"
        if not os.path.exists(pklFile):
            self.printMess("Writing pickle file %s" % pklFile)
            pickle.dump(self, open(pklFile,"w"))
        elif self.force:
            self.printMess("Overwriting pickle file %s" % pklFile)
            pickle.dump(self, open(pklFile,"w"))
        else:
            self.printMess("Pickle file %s already present... doing nothing" % pklFile)

    def getFlagData(self, flag):
        """
            For a given acFileTop flag, return a list of the data related
        """
        block = False
        tFlag = '%FLAG ' + flag
        data = ''

        for rawLine in self.topFileData:
            line = rawLine[:-1]
            if tFlag in line:
                block = True
                continue
            if block and '%FLAG ' in line: break
            if block:
                if '%FORMAT' in line:
                    line = line.strip().strip('%FORMAT()').split('.')[0]
                    for c in line:
                        if c.isalpha():
                            f = int(line.split(c)[1])
                            break
                    continue
                data += line
        # data need format
        sdata = [data[i:i+f].strip() for i in xrange(0,len(data),f)]
        if '+' and '.' in data: # it's a float
            ndata = map(float, sdata)
        elif flag != 'RESIDUE_LABEL':
            try: # try if it's integer
                ndata = map(int, sdata)
            except: # it's string
                ndata = sdata
        else:
            ndata = sdata
        return ndata # a list

    def getResidueLabel(self):
        """
            Get a 3 capital letters code from acFileTop
            Returns a list.
        """
        residueLabel = self.getFlagData('RESIDUE_LABEL')
        residueLabel = map(str, residueLabel)
        if residueLabel[0] != residueLabel[0].upper():
            self.printWarn("residue label '%s' in '%s' is not all UPPERCASE" %
                           (residueLabel[0], self.inputFile))
            self.printWarn("this may raise problem with some applications like CNS")
        self.residueLabel = residueLabel

    def getCoords(self):
        """
            For a given acFileXyz file, return a list of coords as:
            [[x1,y1,z1],[x2,y2,z2], etc.]
        """
        data = ''
        for rawLine in self.xyzFileData[2:]:
            line = rawLine[:-1]
            data += line
        l = len(data)
        ndata = map(float,[data[i:i+12] for i in xrange(0,l,12)])

        gdata = []
        for i in xrange(0, len(ndata), 3):
            gdata.append([ndata[i], ndata[i+1], ndata[i+2]])

        self.printDebug("getCoords done")

        return gdata

    def getAtoms(self):
        """
            Set a list with all atoms objects build from dat in acFileTop
            Set also if molTopol atom type system is gaff or amber
            Set also list atomTypes
            Set also resid
            Set also molTopol total charge
        """
        atomNameList = self.getFlagData('ATOM_NAME')
        atomTypeNameList = self.getFlagData('AMBER_ATOM_TYPE')
        self._atomTypeNameList = atomTypeNameList
        massList = self.getFlagData('MASS')
        chargeList = self.getFlagData('CHARGE')
        resIds = self.getFlagData('RESIDUE_POINTER') + [0]
        #uniqAtomTypeId = self.getFlagData('ATOM_TYPE_INDEX') # for LJ
        balanceChargeList = self.balanceCharges(chargeList)
        coords = self.getCoords()
        ACOEFs, BCOEFs = self.getABCOEFs()

        atoms = []
        atomTypes = []
        tmpList = [] # a list with unique atom types
        totalCharge = 0.0
        countRes = 0
        id = 0
        for atomName in atomNameList:
            atomTypeName = atomTypeNameList[id]
            if id + 1 == resIds[countRes]:
                resid = countRes #self.residueLabel[countRes]
                countRes += 1
            mass = massList[id]
            charge = balanceChargeList[id]
            totalCharge += charge
            coord = coords[id]
            ACOEF = ACOEFs[id]
            BCOEF = BCOEFs[id]
            atomType = AtomType(atomTypeName, mass, ACOEF, BCOEF)
            if atomTypeName not in tmpList:
                tmpList.append(atomTypeName)
                atomTypes.append(atomType)
            atom = Atom(atomName, atomType, id + 1, resid, mass, charge, coord)
            atoms.append(atom)
            id += 1

        if atomTypeName[0].islower():
            self.atomTypeSystem = 'gaff'
        else:
            self.atomTypeSystem = 'amber'

        self.printDebug('Balanced TotalCharge %13.10f' % float(totalCharge / qConv))
        self.totalCharge = int(totalCharge)

        self.atoms = atoms
        self.atomTypes = atomTypes

        self.pbc = None
        if len(coords) == len(atoms) + 2:
            self.pbc = [coords[-2], coords[-1]]
        self.printDebug("PBC = '%s" % self.pbc)
        self.printDebug("getAtoms done")

    def getBonds(self):
        uniqKbList = self.getFlagData('BOND_FORCE_CONSTANT')
        uniqReqList = self.getFlagData('BOND_EQUIL_VALUE')
        bondCodeHList = self.getFlagData('BONDS_INC_HYDROGEN')
        bondCodeNonHList = self.getFlagData('BONDS_WITHOUT_HYDROGEN')
        bondCodeList = bondCodeHList + bondCodeNonHList
        bonds = []
        for i in xrange(0, len(bondCodeList), 3):
            idAtom1 = bondCodeList[i] / 3 # remember python starts with id 0
            idAtom2 = bondCodeList[i+1] / 3
            bondTypeId = bondCodeList[i+2] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            kb = uniqKbList[bondTypeId]
            req = uniqReqList[bondTypeId]
            atoms = [atom1, atom2]
            bond = Bond(atoms, kb, req)
            bonds.append(bond)
        self.bonds = bonds
        self.printDebug("getBonds done")

    def getAngles(self):
        uniqKtList = self.getFlagData('ANGLE_FORCE_CONSTANT')
        uniqTeqList = self.getFlagData('ANGLE_EQUIL_VALUE')
        # for list below, true atom number = index/3 + 1
        angleCodeHList = self.getFlagData('ANGLES_INC_HYDROGEN')
        angleCodeNonHList = self.getFlagData('ANGLES_WITHOUT_HYDROGEN')
        angleCodeList = angleCodeHList + angleCodeNonHList
        angles = []
        for i in xrange(0, len(angleCodeList), 4):
            idAtom1 = angleCodeList[i] / 3 # remember python starts with id 0
            idAtom2 = angleCodeList[i+1] / 3
            idAtom3 = angleCodeList[i+2] / 3
            angleTypeId = angleCodeList[i+3] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            kt = uniqKtList[angleTypeId]
            teq = uniqTeqList[angleTypeId] # angle given in rad in prmtop
            atoms = [atom1, atom2, atom3]
            angle = Angle(atoms, kt, teq)
            angles.append(angle)
        self.angles = angles
        self.printDebug("getAngles done")

    def getDihedrals(self):
        """
            Get dihedrals (proper and imp), condensed list of prop dih and
            atomPairs
        """
        uniqKpList = self.getFlagData('DIHEDRAL_FORCE_CONSTANT')
        uniqPeriodList = self.getFlagData('DIHEDRAL_PERIODICITY')
        uniqPhaseList = self.getFlagData('DIHEDRAL_PHASE')
        # for list below, true atom number = abs(index)/3 + 1
        dihCodeHList = self.getFlagData('DIHEDRALS_INC_HYDROGEN')
        dihCodeNonHList = self.getFlagData('DIHEDRALS_WITHOUT_HYDROGEN')
        dihCodeList = dihCodeHList + dihCodeNonHList
        properDih = []
        improperDih = []
        condProperDih = [] # list of dihedrals condensed by the same quartet
        #atomPairs = []
        atomPairs = set()
        for i in xrange(0, len(dihCodeList), 5):
            idAtom1 = dihCodeList[i] / 3 # remember python starts with id 0
            idAtom2 = dihCodeList[i+1] / 3
            # 3 and 4 indexes can be negative: if id3 < 0, end group interations
            # in amber are to be ignored; if id4 < 0, dihedral is improper
            idAtom3raw = dihCodeList[i+2] / 3 # can be negative -> exclude from 1-4vdw
            idAtom4raw = dihCodeList[i+3] / 3 # can be negative -> Improper
            idAtom3 = abs(idAtom3raw)
            idAtom4 = abs(idAtom4raw)
            dihTypeId = dihCodeList[i+4] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            atom4 = self.atoms[idAtom4]
            kPhi = uniqKpList[dihTypeId] # already divided by IDIVF
            period = int(uniqPeriodList[dihTypeId]) # integer
            phase = uniqPhaseList[dihTypeId]# angle given in rad in prmtop
            atoms = [atom1, atom2, atom3, atom4]
            dihedral = Dihedral(atoms, kPhi, period, phase)
            if idAtom4raw > 0:
                try: atomsPrev = properDih[-1].atoms
                except: atomsPrev = []
                properDih.append(dihedral)
                if idAtom3raw < 0 and atomsPrev == atoms:
                    condProperDih[-1].append(dihedral)
                else:
                    condProperDih.append([dihedral])
                pair = (atom1, atom4)
                #if atomPairs.count(pair) == 0 and idAtom3raw > 0:
                if idAtom3raw > 0:
                    atomPairs.add(pair)
            else:
                improperDih.append(dihedral)
        try: atomPairs = sorted(atomPairs)
        except: pass
        self.properDihedrals = properDih
        self.improperDihedrals = improperDih
        self.condensedProperDihedrals = condProperDih # [[],[],...]
        self.atomPairs = atomPairs # set((atom1, atom2), ...)
        self.printDebug("getDihedrals done")

    def setAtomPairs(self):
        """
            Set a list of pair of atoms pertinent to interaction 1-4 for vdw.
            WRONG: Deprecated
        """
        atomPairs = []
        for item in self.condensedProperDihedrals:
            dih = item[0]
            atom1 = dih.atoms[0]
            atom2 = dih.atoms[3]
            pair = [atom1, atom2]
            if atomPairs.count(pair) == 0:
                atomPairs.append(pair)
        self.atomPairs = atomPairs # [[atom1, atom2], ...]
        self.printDebug("atomPairs done")

    def getExcludedAtoms(self):
        """
            Returns a list of atoms with a list of its excluded atoms up to 3rd
            neighbour.
            It's implicitly indexed, i.e., a sequence of atoms in position n in
            the excludedAtomsList corresponds to atom n (self.atoms) and so on.
            NOT USED
        """
        excludedAtomsIdList = self.getFlagData('EXCLUDED_ATOMS_LIST')
        numberExcludedAtoms = self.getFlagData('NUMBER_EXCLUDED_ATOMS')
        atoms = self.atoms
        interval = 0
        excludedAtomsList = []
        for number in numberExcludedAtoms:
            temp = excludedAtomsIdList[interval:interval + number]
            if temp == [0]:
                excludedAtomsList.append([])
            else:
                excludedAtomsList.append([atoms[a-1] for a in temp])
            interval += number
        self.excludedAtoms = excludedAtomsList
        self.printDebug("getExcludedAtoms")

    def balanceCharges(self, chargeList):
        """
            Note that python is very annoying about floating points.
            Even after balance, there will always be some residue of order e-12
            to e-16, which is believed to vanished once one writes a topology
            file, say, for CNS or GMX, where floats are represented with 4 or 5
            maximum decimals.
        """
        #self.printDebug(chargeList)
        total = sum(chargeList)/qConv
        self.printDebug('charge to be balanced: total %13.10f' % total)
        maxVal = max(chargeList)
        minVal = min(chargeList)
        #print maxVal, minVal
        if abs(maxVal) >= abs(minVal): lim = maxVal
        else: lim = minVal
        limId = chargeList.index(lim)
        diff = (total - round(total)) * qConv
        #print diff
        fix = lim - diff
        chargeList[limId] = fix
        self.printDebug("balanceCharges done")
        return chargeList

    def getABCOEFs(self):
        uniqAtomTypeIdList = self.getFlagData('ATOM_TYPE_INDEX')
        nonBonIdList = self.getFlagData('NONBONDED_PARM_INDEX')
        rawACOEFs = self.getFlagData('LENNARD_JONES_ACOEF')
        rawBCOEFs = self.getFlagData('LENNARD_JONES_BCOEF')
        #print nonBonIdList, len(nonBonIdList), rawACOEFs, len(rawACOEFs)
        ACOEFs = []
        BCOEFs = []
        ntypes = max(uniqAtomTypeIdList)
        #id = 0
        #for atName in self._atomTypeNameList:
        for id in xrange(len(self._atomTypeNameList)):
            #id = self._atomTypeNameList.index(atName)
            atomTypeId = uniqAtomTypeIdList[id]
            index = ntypes * (atomTypeId - 1) + atomTypeId
            nonBondId = nonBonIdList[index - 1]
            #print "*****", index, ntypes, atName, id, atomTypeId, nonBondId
            ACOEFs.append(rawACOEFs[nonBondId - 1])
            BCOEFs.append(rawBCOEFs[nonBondId - 1])
            #id += 1
        #print ACOEFs
        self.printDebug("getABCOEFs done")
        return ACOEFs, BCOEFs

    def setProperDihedralsCoefRB(self):
        """
            It takes self.condensedProperDihedrals and returns
            self.properDihedralsCoefRB, a reduced list of quartet atoms + RB.
            Coeficients ready for GMX (multiplied by 4.184)

            self.properDihedralsCoefRB = [ [atom1,..., atom4], C[0:5] ]

            For proper dihedrals: a quartet of atoms may appear with more than
            one set of parameters and to convert to GMX they are treated as RBs.

            The resulting coefs calculated here may look slighted different from
            the ones calculated by amb2gmx.pl because python is taken full float
            number from prmtop and not rounded numbers from rdparm.out as
            amb2gmx.pl does.
        """
        properDihedralsCoefRB = []
        for item in self.condensedProperDihedrals:
            V = 6 * [0.0]
            C = 6 * [0.0]
            for dih in item:
                period = dih.period
                kPhi = dih.kPhi # in rad
                phase = int(dih.phase * radPi)
                if kPhi > 0: V[period] = 2 * kPhi * cal
                if period == 1:
                    C[0] += 0.5 * V[period]
                    if phase == 0:
                        C[1] -= 0.5 * V[period]
                    else:
                        C[1] += 0.5 * V[period]
                elif period == 2:
                    if phase == 180:
                        C[0] += V[period]
                        C[2] -= V[period]
                    else:
                        C[2] += V[period]
                elif period == 3:
                    C[0] += 0.5 * V[period]
                    if phase == 0:
                        C[1] += 1.5 * V[period]
                        C[3] -= 2 * V[period]
                    else:
                        C[1] -= 1.5 * V[period]
                        C[3] += 2 * V[period]
                elif period == 4:
                    if phase == 180:
                        C[2] += 4 * V[period]
                        C[4] -= 4 * V[period]
                    else:
                        C[0] += V[period]
                        C[2] -= 4 * V[period]
                        C[4] += 4 * V[period]
                #print kPhi, period, phase, V, C
            properDihedralsCoefRB.append([item[0].atoms,C])

        self.printDebug("setProperDihedralsCoefRB done")

        self.properDihedralsCoefRB = properDihedralsCoefRB

    def writeCharmmTopolFiles(self):

        self.printMess("Writing CHARMM files\n")

        #self.makeDir()

        at = self.atomType
        self.getResidueLabel()
        res = self.resName #self.residueLabel[0]
        #print res, self.residueLabel, type(self.residueLabel)

        cmd = '%s -i %s -fi mol2 -o %s -fo charmm -s 2 -df 0 -at %s \
        -pf y -rn %s' % (self.acExe, self.acMol2FileName, self.charmmBase, at, res)

        if self.debug:
            cmd = cmd.replace('-pf y', '-pf n')
        _log = getoutput(cmd)

    def writePdb(self, file):
        """
            Write a new PDB file with the atom names defined by Antechamber
            Input: file path string
            The format generated here use is slightly different from
            http://www.wwpdb.org/documentation/format23/sect9.html respected to
            atom name
        """
        #TODO: assuming only one residue ('1')
        pdbFile = open(file, 'w')
        fbase = os.path.basename(file)
        pdbFile.write("REMARK "+ head % (fbase, date))
        id = 1
        for atom in self.atoms:
            #id = self.atoms.index(atom) + 1
            aName = atom.atomName
            if len(aName) != 4:
                aName = ' ' + aName
            s = aName[1]
            rName = self.residueLabel[0]
            x = atom.coords[0]
            y = atom.coords[1]
            z = atom.coords[2]
            line = "%-6s%5d %-5s%3s Z%4d%s%8.3f%8.3f%8.3f%6.2f%6.2f%s%2s\n" % \
            ('ATOM', id, aName, rName, 1, 4*' ', x, y, z, 1.0, 0.0, 10*' ', s)
            pdbFile.write(line)
            id += 1
        pdbFile.write('END\n')

    def writeGromacsTopolFiles(self, amb2gmx = False):
        """
            # from ~/Programmes/amber10/dat/leap/parm/gaff.dat
            #atom type        atomic mass        atomic polarizability        comments
            ca                12.01                 0.360                    Sp2 C in pure aromatic systems
            ha                1.008                 0.135                    H bonded to aromatic carbon

            #bonded atoms        harmonic force kcal/mol/A^2       eq. dist. Ang.  comments
            ca-ha                  344.3*                           1.087**         SOURCE3  1496    0.0024    0.0045
            * for gmx: 344.3 * 4.184 * 100 * 2 = 288110 kJ/mol/nm^2 (why factor 2?)
            ** convert Ang to nm ( div by 10) for gmx: 1.087 A = 0.1087 nm
            # CA HA         1    0.10800   307105.6 ; ged from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from ffamber99bon.itp)
            # CA-HA  367.0    1.080       changed from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from parm99.dat)

            # angle        HF kcal/mol/rad^2    eq angle degrees     comments
            ca-ca-ha        48.5*             120.01                SOURCE3 2980   0.1509   0.2511
            * to convert to gmx: 48.5 * 4.184 * 2 = 405.848 kJ/mol/rad^2 (why factor 2?)
            # CA  CA  HA           1   120.000    418.400 ; new99 (from ffamber99bon.itp)
            # CA-CA-HA    50.0      120.00 (from parm99.dat)

            # dihedral    idivf        barrier hight/2 kcal/mol  phase degrees       periodicity     comments
            X -ca-ca-X    4           14.500*                     180.000                2.000             intrpol.bsd.on C6H6
            * to convert to gmx: 14.5/4 * 4.184 * 2 (?) (yes in amb2gmx, no in topolbuild, why?) = 30.334 or 15.167 kJ/mol
            # X -CA-CA-X    4   14.50        180.0             2.         intrpol.bsd.on C6H6 (from parm99.dat)
            # X   CA  CA  X     3    30.33400     0.00000   -30.33400     0.00000     0.00000     0.00000   ; intrpol.bsd.on C6H6
            ;propers treated as RBs in GROMACS to use combine multiple AMBER torsions per quartet (from ffamber99bon.itp)

            # impr. dihedral        barrier hight/2      phase degrees       periodicity     comments
            X -X -ca-ha             1.1*                  180.                      2.                   bsd.on C6H6 nmodes
            * to convert to gmx: 1.1 * 4.184 = 4.6024 kJ/mol/rad^2
            # X -X -CA-HA         1.1          180.          2.           bsd.on C6H6 nmodes (from parm99.dat)
            # X   X   CA  HA       1      180.00     4.60240     2      ; bsd.on C6H6 nmodes
            ;impropers treated as propers in GROMACS to use correct AMBER analytical function (from ffamber99bon.itp)

            # 6-12 parms     sigma = 2 * r * 2^(-1/6)    epsilon
            # atomtype        radius Ang.                    pot. well depth kcal/mol      comments
              ha                  1.4590*                      0.0150**                         Spellmeyer
              ca                  1.9080                    0.0860                            OPLS
            * to convert to gmx:
                sigma = 1.4590 * 2^(-1/6) * 2 = 2 * 1.29982 Ang. = 2 * 0.129982 nm  = 1.4590 * 2^(5/6)/10 =  0.259964 nm
            ** to convert to gmx: 0.0150 * 4.184 = 0.06276 kJ/mol
            # amber99_3    CA     0.0000  0.0000  A   3.39967e-01  3.59824e-01 (from ffamber99nb.itp)
            # amber99_22   HA     0.0000  0.0000  A   2.59964e-01  6.27600e-02 (from ffamber99nb.itp)
            # C*          1.9080  0.0860             Spellmeyer
            # HA          1.4590  0.0150             Spellmeyer (from parm99.dat)
            # to convert r and epsilon to ACOEF and BCOEF
            # ACOEF = sqrt(e1*e2) * (r1 + r2)^12 ; BCOEF = 2 * sqrt(e1*e2) * (r1 + r2)^6 = 2 * ACOEF/(r1+r2)^6
            # to convert ACOEF and BCOEF to r and espsilon
            # r = 0.5 * (2*ACOEF/BCOEF)^(1/6); ep = BCOEF^2/(4*ACOEF)
            # to convert ACOEF and BCOEF to sigma and epsilon (GMX)
            # sigma = (ACOEF/BCOEF)^(1/6) * 0.1 ; epsilon = 4.184 * BCOEF^2/(4*ACOEF)
            #   ca   ca       819971.66        531.10
            #   ca   ha        76245.15        104.66
            #   ha   ha         5716.30         18.52

            For proper dihedrals: a quartet of atoms may appear with more than
            one set of parameters and to convert to GMX they are treated as RBs;
            use the algorithm:
              for(my $j=$i;$j<=$lines;$j++){
                my $period = $pn{$j};
                if($pk{$j}>0) {
                  $V[$period] = 2*$pk{$j}*$cal;
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
        """

        self.printMess("Writing GROMACS files\n")

        self.setAtomType4Gromacs()

        self.writeGroFile()

        self.writeGromacsTop(amb2gmx = amb2gmx)

        self.writeMdpFiles()

    def setAtomType4Gromacs(self):
        """Atom types names in Gromacs TOP file are not case sensitive;
           this routine will append a '_' to lower case atom type.
           E.g.: CA and ca -> CA and ca_
        """
        atNames = [at.atomTypeName for at in self.atomTypes]
        #print atNames
        delAtomTypes = []
        modAtomTypes = []
        atomTypesGromacs = []
        dictAtomTypes = {}
        for at in self.atomTypes:
            atName = at.atomTypeName
            dictAtomTypes[atName] = at
            if atName.islower() and atName.upper() in atNames:
                #print atName, atName.upper()
                atUpper = self.atomTypes[atNames.index(atName.upper())]
                #print at.atomTypeName,at.mass, at.ACOEF, at.BCOEF
                #print atUpper.atomTypeName, atUpper.mass, atUpper.ACOEF, atUpper.BCOEF
                if at.ACOEF is atUpper.ACOEF and at.BCOEF is at.BCOEF:
                    delAtomTypes.append(atName)
                else:
                    newAtName = atName+'_'
                    modAtomTypes.append(atName)
                    atomType = AtomType(newAtName, at.mass, at.ACOEF, at.BCOEF)
                    atomTypesGromacs.append(atomType)
                    dictAtomTypes[newAtName] = atomType
            else:
                atomTypesGromacs.append(at)

        atomsGromacs = []
        for a in self.atoms:
            atName = a.atomType.atomTypeName
            if atName in delAtomTypes:
                atom = Atom(a.atomName, dictAtomTypes[atName.upper()], a.id, \
                            a.resid, a.mass, a.charge, a.coords)
                atomsGromacs.append(atom)
            elif atName in modAtomTypes:
                atom = Atom(a.atomName, dictAtomTypes[atName + '_'], a.id, \
                            a.resid, a.mass, a.charge, a.coords)
                atomsGromacs.append(atom)
            else:
                atomsGromacs.append(a)

        self.atomTypesGromacs = atomTypesGromacs
        self.atomsGromacs = atomsGromacs
        #print [i.atomTypeName for i in atomTypesGromacs]
        #print modAtomTypes
        #print delAtomTypes

    def writeGromacsTop(self, amb2gmx = False):
        if self.atomTypeSystem == 'amber':
            d2opls = dictAtomTypeAmb2OplsGmxCode
        else:
            d2opls = dictAtomTypeGaff2OplsGmxCode

        topText = []
        itpText = []
        oitpText = []
        otopText = []
        top = self.baseName+'_GMX.top'
        itp = self.baseName+'_GMX.itp'
        otop = self.baseName+'_GMX_OPLS.top'
        oitp = self.baseName+'_GMX_OPLS.itp'

        headDefault = \
"""
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333
"""
        headItp = \
"""
; Include %s topology
#include "%s"
"""
        headOpls = \
"""
; Include forcefield parameters
#include "ffoplsaa.itp"
"""
        headSystem = \
"""
[ system ]
 %s
"""
        headMols = \
"""
[ molecules ]
; Compound        nmols
"""
        headAtomtypes = \
"""
[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
"""
        headAtomtypesOpls = \
"""
; For OPLS atomtypes manual fine tuning
; AC_at:OPLS_at:OPLScode: Possible_Aternatives (see ffoplsaa.atp and ffoplsaanb.itp)
"""
        headMoleculetype = \
"""
[ moleculetype ]
;name            nrexcl
 %-16s 3
"""
        headAtoms = \
"""
[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
"""
        headBonds = \
"""
[ bonds ]
;   ai     aj funct   r             k
"""
        headPairs = \
"""
[ pairs ]
;   ai     aj    funct
"""
        headAngles = \
"""
[ angles ]
;   ai     aj     ak    funct   theta         cth
"""
        headProDih = \
"""
[ dihedrals ] ; propers
; treated as RBs in GROMACS to use combine multiple AMBER torsions per quartet
; i   j   k   l func   C0        C1        C2        C3        C4        C5
"""
        headImpDih = \
"""
[ dihedrals ] ; impropers
; treated as propers in GROMACS to use correct AMBER analytical function
; i   j   k   l func  phase     kd      pn
"""
        _headTopWaterTip3p = \
"""
[ bondtypes ]
  ; i    j      func       b0          kb
  OW    HW         1    0.09572   462750.4 ; TIP3P water
  HW    HW         1    0.15139   462750.4 ; TIP3P water

[ angletypes ]
  ;  i    j    k  func       th0       cth
  HW  OW  HW           1   104.520    836.800 ; TIP3P water
  HW  HW  OW           1   127.740      0.000 ; (found in crystallographic water with 3 bonds)
"""

        _headTopWaterSpce = \
"""
[ bondtypes ]
  ; i    j      func       b0          kb
  OW    HW         1    0.1       462750.4 ; SPCE water
  HW    HW         1    0.1633    462750.4 ; SPCE water

[ angletypes ]
  ;  i    j    k  func       th0       cth
  HW  OW  HW           1   109.47      836.800 ; SPCE water
  HW  HW  OW           1   125.265     0.000 ; SPCE water
"""
#NOTE: headTopWaterTip3p and headTopWaterSpce actually do NOTHING
        headCl = \
"""
[ moleculetype ]
  ; molname       nrexcl
  Na+             1

[ atoms ]
  ; id    at type res nr  residu name     at name  cg nr  charge   mass
    1       IP      1          Na+         Na+       1      1     22.9898
"""
        headNa = \
"""
[ moleculetype ]
  ; molname       nrexcl
  Cl-             1

[ atoms ]
  ; id    at type res nr  residu name     at name  cg nr  charge   mass
    1       IM      1         Cl-           Cl-      1     -1     35.45300
"""

        headWaterTip3p = \
"""
[ moleculetype ]
; molname       nrexcl ; TIP3P model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.834   16.00000
     2     HW      1     WAT    H1      1      0.417    1.00800
     3     HW      1     WAT    H2      1      0.417    1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.09572   462750.4 0.09572   462750.4
1   3   1   0.09572   462750.4 0.09572   462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   104.520    836.800  104.520    836.800
#else
[ settles ]
; i j   funct   length
1   1   0.09572 0.15139

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
"""

        headWaterSpce = \
"""
[ moleculetype ]
; molname       nrexcl ; SPCE model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.8476  15.99940
     2     HW      1     WAT    H1      1      0.4238   1.00800
     3     HW      1     WAT    H2      1      0.4238   1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.1 462750.4  0.1     462750.4
1   3   1   0.1 462750.4  0.1     462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   109.47  836.800 109.47  836.800
#else
[ settles ]
; OW    funct   doh dhh
1   1   0.1 0.16330

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
"""
#NOTE: headWaterTip3p and headWaterSpce actually do the real thing
#      so, skipping headTopWaterTip3p and headWaterTip3p
        #headTopWater = headTopWaterTip3p
        headWater = headWaterTip3p

        nWat = nCl = nNa = 0
        #topFile.write("; " + head % (top, date))
        topText.append("; " + head % (top, date))
        otopText.append("; " + head % (otop, date))
        #topFile.write(headDefault)
        topText.append(headDefault)

        nSolute = 0
        if not amb2gmx:
            topText.append(headItp % (itp, itp))
            otopText.append(headOpls)
            otopText.append(headItp % (itp, itp))
            itpText.append("; " + head % (itp, date))
            oitpText.append("; " + head % (oitp, date))

        self.printDebug("atomTypes %i" % len(self.atomTypesGromacs))
        temp = []
        otemp = []
        for aType in self.atomTypesGromacs:
            aTypeName = aType.atomTypeName
            oaCode = d2opls.get(aTypeName,['x','0'])[:-1]
            aTypeNameOpls = oplsCode2AtomTypeDict.get(oaCode[0],'x')
            A = aType.ACOEF
            B = aType.BCOEF
            # one cannot infer sigma or epsilon for B = 0, assuming 0 for them
            if B == 0.0:
                sigma, epsilon, r0, epAmber = 0, 0, 0, 0
            else:
                r0 = 0.5 * math.pow((2*A/B), (1.0/6))
                epAmber = 0.25 * B*B/A
                sigma = 0.1 * math.pow((A/B), (1.0/6))
                epsilon = cal * epAmber
            if aTypeName == 'OW':
                if A == 629362.166 and B == 625.267765:
                    #headTopWater = headTopWaterSpce
                    headWater = headWaterSpce
            # OW 629362.166 625.267765 spce
            # OW 581935.564 594.825035 tip3p
            #       print aTypeName, A, B
            line = " %-8s %-11s %3.5f  %3.5f   A   %13.5e %13.5e" % \
            (aTypeName, aTypeName, 0.0, 0.0, sigma, epsilon) + \
            " ; %4.2f  %1.4f\n" % (r0, epAmber)
            oline = "; %s:%s:opls_%s: %s\n" % (aTypeName,aTypeNameOpls,oaCode[0],`oaCode[1:]`)
            #tmpFile.write(line)
            temp.append(line)
            otemp.append(oline)
        if amb2gmx:
            topText.append(headAtomtypes)
            topText += temp
            nWat = self.residueLabel.count('WAT')
            nCl = self.residueLabel.count('Cl-')
            nNa = self.residueLabel.count('Na+')
        else:
            itpText.append(headAtomtypes)
            itpText += temp
            oitpText.append(headAtomtypesOpls)
            oitpText += otemp
        self.printDebug("GMX atomtypes done")

        if len(self.atoms) > 3 * nWat + nCl + nNa:
            nSolute = 1

        if nWat:
            #topText.append(headTopWater)
            self.printDebug("type of water '%s'" % headWater[43:48].strip())

        if nSolute:
            if amb2gmx:
                topText.append(headMoleculetype % self.baseName)
            else:
                itpText.append(headMoleculetype % self.baseName)
                oitpText.append(headMoleculetype % self.baseName)

        self.printDebug("atoms %i" % len(self.atoms))
        qtot = 0.0
        count = 1
        temp = []
        otemp = []
        id2oplsATDict = {}
        for atom in self.atomsGromacs:
            resid = atom.resid
            resname = self.residueLabel[resid]
            if resname in ['Cl-', 'Na+', 'WAT' ]:
                break
            aName = atom.atomName
            aType = atom.atomType.atomTypeName
            oItem = d2opls.get(aType,['x',0])
            oplsAtName = oplsCode2AtomTypeDict.get(oItem[0],'x')
            id = atom.id
            id2oplsATDict[id] = oplsAtName
            oaCode = 'opls_'+oItem[0]
            charge = atom.charge
            mass = atom.mass
            omass = float(oItem[-1])
            qtot += charge
            resnr = resid + 1
            line = "%6d %4s %5d %5s %5s %4d %12.5f %12.5f ; qtot %1.3f\n" % \
            (count, aType, resnr, resname, aName, count, charge, mass, qtot)
            oline = "%6d %4s %5d %5s %5s %4d %12.5f %12.5f ; qtot % 3.3f  %-4s\n" % \
            (count, oaCode, resnr, resname, aName, count, charge, omass, qtot, oplsAtName)
            count += 1
            temp.append(line)
            otemp.append(oline)
        if temp:
            if amb2gmx:
                topText.append(headAtoms)
                topText += temp
            else:
                itpText.append(headAtoms)
                itpText += temp
                oitpText.append(headAtoms)
                oitpText += otemp
        self.printDebug("GMX atoms done")

        # remove bond of water
        self.printDebug("bonds %i" % len(self.bonds))
        temp = []
        otemp = []
        for bond in self.bonds:
            res1 = self.residueLabel[bond.atoms[0].resid]
            res2 = self.residueLabel[bond.atoms[0].resid]
            if 'WAT' in [res1, res2] : continue
            a1Name = bond.atoms[0].atomName
            a2Name = bond.atoms[1].atomName
            id1 = bond.atoms[0].id
            id2 = bond.atoms[1].id
            oat1 = id2oplsATDict[id1]
            oat2 = id2oplsATDict[id2]
            line = "%6i %6i %3i %13.4e %13.4e ; %6s - %-6s\n" % (id1, id2, 1,
                   bond.rEq * 0.1, bond.kBond * 200 * cal, a1Name, a2Name)
            oline = "%6i %6i %3i ; %13.4e %13.4e ; %6s - %-6s %6s - %-6s\n" % \
                (id1, id2, 1, bond.rEq * 0.1, bond.kBond * 200 * cal, a1Name,
                 a2Name, oat1, oat2)
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if amb2gmx:
                topText.append(headBonds)
                topText += temp
            else:
                itpText.append(headBonds)
                itpText += temp
                oitpText.append(headBonds)
                oitpText += otemp
        self.printDebug("GMX bonds done")

        self.printDebug("atomPairs %i" % len(self.atomPairs))
        temp = []
        for pair in self.atomPairs:
            #if not printed:
            #    tmpFile.write(headPairs)
            #    printed = True
            a1Name = pair[0].atomName
            a2Name = pair[1].atomName
            id1 = pair[0].id
            id2 = pair[1].id
            #id1 = self.atoms.index(pair[0]) + 1
            #id2 = self.atoms.index(pair[1]) + 1
            line = "%6i %6i %6i ; %6s - %-6s\n" % (id1, id2, 1, a1Name,
                                                        a2Name)
            temp.append(line)
        temp.sort()
        if temp:
            if amb2gmx:
                topText.append(headPairs)
                topText += temp
            else:
                itpText.append(headPairs)
                itpText += temp
                oitpText.append(headPairs)
                oitpText += temp
        self.printDebug("GMX pairs done")

        self.printDebug("angles %i" % len(self.angles))
        temp = []
        otemp = []
        for angle in self.angles:
            a1 = angle.atoms[0].atomName
            a2 = angle.atoms[1].atomName
            a3 = angle.atoms[2].atomName
            id1 = angle.atoms[0].id
            id2 = angle.atoms[1].id
            id3 = angle.atoms[2].id
            oat1 = id2oplsATDict[id1]
            oat2 = id2oplsATDict[id2]
            oat3 = id2oplsATDict[id3]
            line = "%6i %6i %6i %6i %13.4e %13.4e ; %6s - %-6s - %-6s\n" % (id1, id2,
            id3, 1, angle.thetaEq * radPi, 2 * cal * angle.kTheta, a1, a2, a3)
            oline = "%6i %6i %6i %6i ; %13.4e %13.4e ; %6s - %-4s - %-6s %4s - %+4s - %-4s\n" % \
                (id1, id2, id3, 1, angle.thetaEq * radPi, 2 * cal * angle.kTheta,
                 a1, a2, a3, oat1, oat2, oat3)
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if amb2gmx:
                topText.append(headAngles)
                topText += temp
            else:
                itpText.append(headAngles)
                itpText += temp
                oitpText.append(headAngles)
                oitpText += otemp
        self.printDebug("GMX angles done")

        self.setProperDihedralsCoefRB()
        self.printDebug("properDihedralsCoefRB %i" % len(self.properDihedralsCoefRB))
        temp = []
        otemp = []
        for dih in self.properDihedralsCoefRB:
            a1 = dih[0][0].atomName
            a2 = dih[0][1].atomName
            a3 = dih[0][2].atomName
            a4 = dih[0][3].atomName
            #id1 = self.atoms.index(dih[0][0]) + 1
            id1 = dih[0][0].id
            #id2 = self.atoms.index(dih[0][1]) + 1
            id2 = dih[0][1].id
            #id3 = self.atoms.index(dih[0][2]) + 1
            id3 = dih[0][2].id
            #id4 = self.atoms.index(dih[0][3]) + 1
            id4 = dih[0][3].id
            oat1 = id2oplsATDict[id1]
            oat2 = id2oplsATDict[id2]
            oat3 = id2oplsATDict[id3]
            oat4 = id2oplsATDict[id4]
            c0, c1, c2, c3, c4, c5 = dih[1]
            line = \
            "%3i %3i %3i %3i %3i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" %\
            (id1, id2, id3, id4, 3, c0, c1, c2, c3, c4, c5) \
            + " ; %6s-%6s-%6s-%6s\n" % (a1, a2, a3, a4)
            oline = \
            "%3i %3i %3i %3i %3i ; %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" %\
            (id1, id2, id3, id4, 3, c0, c1, c2, c3, c4, c5) \
            + " ; %6s-%6s-%6s-%6s    %4s-%4s-%4s-%4s\n" % (a1, a2, a3, a4, oat1, oat2, oat3, oat4)
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if amb2gmx:
                topText.append(headProDih)
                topText += temp
            else:
                itpText.append(headProDih)
                itpText += temp
                oitpText.append(headProDih)
                oitpText += otemp
        self.printDebug("GMX proper dihedrals done")

        self.printDebug("improperDihedrals %i" % len(self.improperDihedrals))
        temp = []
        otemp = []
        for dih in self.improperDihedrals:
            a1 = dih.atoms[0].atomName
            a2 = dih.atoms[1].atomName
            a3 = dih.atoms[2].atomName
            a4 = dih.atoms[3].atomName
            #id1 = self.atoms.index(dih.atoms[0]) + 1
            id1 = dih.atoms[0].id
            #id2 = self.atoms.index(dih.atoms[1]) + 1
            id2 = dih.atoms[1].id
            #id3 = self.atoms.index(dih.atoms[2]) + 1
            id3 = dih.atoms[2].id
            #id4 = self.atoms.index(dih.atoms[3]) + 1
            id4 = dih.atoms[3].id
            kd = dih.kPhi * cal
            pn = dih.period
            ph = dih.phase * radPi
            line = "%3i %3i %3i %3i %3i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % \
                            (id1, id2, id3, id4, 1, ph, kd, pn, a1, a2, a3, a4)
            oline = "%3i %3i %3i %3i %3i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % \
                            (id1, id2, id3, id4, 1, ph, kd, pn, a1, a2, a3, a4)
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if amb2gmx:
                topText.append(headImpDih)
                topText += temp
            else:
                itpText.append(headImpDih)
                itpText += temp
                oitpText.append(headImpDih)
                oitpText += otemp
        self.printDebug("GMX improper dihedrals done")

        if nNa:
            topText.append(headNa)
        if nCl:
            topText.append(headCl)
        if nWat:
            topText.append(headWater)

        topText.append(headSystem % (self.baseName))
        topText.append(headMols)
        otopText.append(headSystem % (self.baseName))
        otopText.append(headMols)

        if nSolute > 0:
            topText.append(" %-16s %-6i\n" % (self.baseName, nSolute))
            otopText.append(" %-16s %-6i\n" % (self.baseName, nSolute))

        if nNa:
            topText.append(" %-16s %-6i\n" % ('Na+', nNa))
        if nCl:
            topText.append(" %-16s %-6i\n" % ('Cl-', nCl))
        if nWat:
            topText.append(" %-16s %-6i\n" % ('WAT', nWat))

        gmxDir = os.path.abspath('.')
        topFileName = os.path.join(gmxDir, top)
        topFile = open(topFileName, 'w')
        topFile.writelines(topText)

        if not amb2gmx:
            itpFileName = os.path.join(gmxDir, itp)
            itpFile = open(itpFileName, 'w')
            itpFile.writelines(itpText)
            oitpFileName = os.path.join(gmxDir, oitp)
            oitpFile = open(oitpFileName, 'w')
            oitpFile.writelines(oitpText)
            otopFileName = os.path.join(gmxDir, otop)
            otopFile = open(otopFileName, 'w')
            otopFile.writelines(otopText)

    def writeGroFile(self):
        #print "Writing GROMACS GRO file\n"
        self.printDebug("writing GRO file")
        gro = self.baseName+'_GMX.gro'
        gmxDir = os.path.abspath('.')
        groFileName = os.path.join(gmxDir, gro)
        groFile = open(groFileName, 'w')
        groFile.write(head % (gro, date))
        groFile.write(" %i\n" % len(self.atoms))
        count = 1
        for atom in self.atoms:
            coords = [c * 0.1 for c in atom.coords]
            resid = atom.resid
            line = "%5d%-4s%6s%5d%8.3f%8.3f%8.3f\n" % \
                   (resid+1, self.residueLabel[resid], atom.atomName,
                    count, coords[0], coords[1], coords[2])
            count += 1
            groFile.write(line)
        if self.pbc:
            boxX = self.pbc[0][0] * 0.1
            boxY = self.pbc[0][1] * 0.1
            boxZ = self.pbc[0][2] * 0.1
            vX = self.pbc[1][0]
            #vY = self.pbc[1][1]
            #vZ = self.pbc[1][2]
            if vX == 90.0:
                self.printDebug("PBC triclinic")
                text = "%11.5f %11.5f %11.5f\n" % (boxX, boxY, boxZ)
            elif round(vX,2) == 109.47:
                self.printDebug("PBC octahedron")
                f1 = 0.471405 # 1/3 * sqrt(2)
                f2 = 0.333333 * boxX
                v22 = boxY * 2 * f1
                v33 = boxZ * f1 * 1.73205 # f1 * sqrt(3)
                v21 = v31 = v32 = 0.0
                v12 = f2
                v13 = - f2
                v23 = f1 * boxX
                text = "%11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n" % \
                   (boxX, v22, v33, v21, v31, v12, v32, v13, v23)
        else:
            self.printDebug("Box size estimated")
            X = [a.coords[0] * 0.1 for a in self.atoms]
            Y = [a.coords[1] * 0.1 for a in self.atoms]
            Z = [a.coords[2] * 0.1 for a in self.atoms]
            boxX = max(X) - min(X) #+ 2.0 # 2.0 is double of rlist
            boxY = max(Y) - min(Y) #+ 2.0
            boxZ = max(Z) - min(Z) #+ 2.0
            text = "%11.5f %11.5f %11.5f\n" % (boxX, boxY, boxZ)
        groFile.write(text)

    def writeMdpFiles(self):
        emMdp = \
"""
cpp                      = /usr/bin/cpp
define                   = -DFLEXIBLE
integrator               = steep
nsteps                   = 500
constraints              = none
emtol                    = 1000.0
emstep                   = 0.01
nstcomm                  = 1
ns_type                  = simple
nstlist                  = 0
rlist                    = 0
rcoulomb                 = 0
rvdw                     = 0
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstxout                  = 1
pbc                      = no
"""
        mdMdp = \
"""
cpp                      = /usr/bin/cpp
define                   = ;-DFLEXIBLE
integrator               = md
nsteps                   = 500
constraints              = none
emtol                    = 1000.0
emstep                   = 0.01
comm_mode                = angular
ns_type                  = simple
nstlist                  = 0
rlist                    = 0
rcoulomb                 = 0
rvdw                     = 0
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstxout                  = 1
pbc                      = no
"""
        emMdpFile = open('em.mdp', 'w')
        mdMdpFile = open('md.mdp', 'w')
        emMdpFile.write(emMdp)
        mdMdpFile.write(mdMdp)

    def writeCnsTopolFiles(self):
        autoAngleFlag = True
        autoDihFlag   = True
        cnsDir = os.path.abspath('.')

        pdb = self.baseName+'_NEW.pdb'
        par = self.baseName+'_CNS.par'
        top = self.baseName+'_CNS.top'
        inp = self.baseName+'_CNS.inp'

        pdbFileName = os.path.join(cnsDir, pdb)
        parFileName = os.path.join(cnsDir, par)
        topFileName = os.path.join(cnsDir, top)
        inpFileName = os.path.join(cnsDir, inp)

        self.CnsTopFileName = topFileName
        self.CnsInpFileName = inpFileName
        self.CnsParFileName = parFileName
        self.CnsPdbFileName = pdbFileName

        parFile = open(parFileName, 'w')
        topFile = open(topFileName, 'w')
        inpFile = open(inpFileName, 'w')

        self.printMess("Writing NEW PDB file\n")
        self.writePdb(pdbFileName)

        self.printMess("Writing CNS/XPLOR files\n")

        #print "Writing CNS PAR file\n"
        parFile.write("Remarks " + head % (par, date))
        parFile.write("\nset echo=false end\n")

        parFile.write("\n{ Bonds: atomType1 atomType2 kb r0 }\n")
        lineSet = []
        for bond in self.bonds:
            a1Type = bond.atoms[0].atomType.atomTypeName+'_'
            a2Type = bond.atoms[1].atomType.atomTypeName+'_'
            kb = 1000.0
            if not self.allhdg:
                kb = bond.kBond
            r0 = bond.rEq
            line = "BOND %5s %5s %8.1f %8.4f\n" % (a1Type, a2Type, kb, r0)
            lineRev = "BOND %5s %5s %8.1f %8.4f\n" % (a2Type, a1Type, kb, r0)
            if line not in lineSet:
                if lineRev not in lineSet:
                    lineSet.append(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Angles: aType1 aType2 aType3 kt t0 }\n")
        lineSet = []
        for angle in self.angles:
            a1 = angle.atoms[0].atomType.atomTypeName+'_'
            a2 = angle.atoms[1].atomType.atomTypeName+'_'
            a3 = angle.atoms[2].atomType.atomTypeName+'_'
            kt = 500.0
            if not self.allhdg:
                kt = angle.kTheta
            t0 = angle.thetaEq * radPi
            line = "ANGLe %5s %5s %5s %8.1f %8.2f\n" % (a1, a2, a3, kt, t0)
            lineRev = "ANGLe %5s %5s %5s %8.1f %8.2f\n" % (a3, a2, a1, kt, t0)
            if line not in lineSet:
                if lineRev not in lineSet:
                    lineSet.append(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Proper Dihedrals: aType1 aType2 aType3 aType4 kt per\
iod phase }\n")
        lineSet = set()
        for item in self.condensedProperDihedrals:
            seq = ''
            id = 0
            for dih in item:
                #id = item.index(dih)
                l = len(item)
                a1 = dih.atoms[0].atomType.atomTypeName+'_'
                a2 = dih.atoms[1].atomType.atomTypeName+'_'
                a3 = dih.atoms[2].atomType.atomTypeName+'_'
                a4 = dih.atoms[3].atomType.atomTypeName+'_'
                kp = 750.0
                if not self.allhdg:
                    kp = dih.kPhi
                p = dih.period
                ph = dih.phase * radPi
                if l > 1:
                    if id == 0:
                        line = "DIHEdral %5s %5s %5s %5s  MULT %1i %7.3f %4i %8\
.2f\n" % (a1, a2, a3, a4, l, kp, p, ph)
                    else:
                        line = "%s %7.3f %4i %8.2f\n" % (40*" ", kp, p, ph)
                else:
                    line = "DIHEdral %5s %5s %5s %5s %15.3f %4i %8.2f\n" % (a1,
                                                          a2, a3, a4, kp, p, ph)
                seq += line
                id += 1
            lineSet.add(seq)
        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Improper Dihedrals: aType1 aType2 aType3 aType4 kt p\
eriod phase }\n")
        lineSet = set()
        for idh in self.improperDihedrals:
            a1 = idh.atoms[0].atomType.atomTypeName+'_'
            a2 = idh.atoms[1].atomType.atomTypeName+'_'
            a3 = idh.atoms[2].atomType.atomTypeName+'_'
            a4 = idh.atoms[3].atomType.atomTypeName+'_'
            kp = 750.0
            if not self.allhdg:
                kp = idh.kPhi
            p = idh.period
            ph = idh.phase * radPi
            line = "IMPRoper %5s %5s %5s %5s %13.1f %4i %8.2f\n" % (a1, a2, a3,
                                                                 a4, kp, p, ph)
            lineSet.add(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Nonbonded: Type Emin sigma; (1-4): Emin/2 sigma }\n")
        for at in self.atomTypes:
            A = at.ACOEF
            B = at.BCOEF
            atName = at.atomTypeName+'_'
            if B == 0.0:
                sigma = epAmber = ep2 = sig2 = 0.0
            else:
                epAmber = 0.25 * B*B/A
                ep2 = epAmber/2.0
                sigma = math.pow((A/B), (1.0/6))
                sig2 = sigma
            line = "NONBonded %5s %11.6f %11.6f %11.6f %11.6f\n" % (atName,
                    epAmber, sigma, ep2, sig2)
            parFile.write(line)
        parFile.write("\nset echo=true end\n")

        #print "Writing CNS TOP file\n"
        topFile.write("Remarks " + head % (top, date))
        topFile.write("\nset echo=false end\n")
        topFile.write("\nautogenerate angles=%s dihedrals=%s end\n" %
                      (autoAngleFlag, autoDihFlag))

        topFile.write("\n{ atomType  mass }\n")
        for at in self.atomTypes:
            atType = at.atomTypeName+'_'
            mass = at.mass
            line = "MASS %-5s %8.3f\n" % (atType, mass)
            topFile.write(line)

        topFile.write("\nRESIdue %s\n" % self.residueLabel[0])
        topFile.write("\nGROUP\n")

        topFile.write("\n{ atomName  atomType  Charge }\n")
        for at in self.atoms:
            atName = at.atomName
            atType = at.atomType.atomTypeName+'_'
            charge = at.charge
            line = "ATOM %-5s TYPE= %-5s CHARGE= %8.4f END\n" % (atName, atType,
                                                                 charge)
            topFile.write(line)

        topFile.write("\n{ Bonds: atomName1  atomName2 }\n")
        for bond in self.bonds:
            a1Name = bond.atoms[0].atomName
            a2Name = bond.atoms[1].atomName
            line = "BOND %-5s %-5s\n" % (a1Name, a2Name)
            topFile.write(line)

        if not autoAngleFlag:
            topFile.write("\n{ Angles: atomName1 atomName2 atomName3}\n")
            for angle in self.angles:
                a1Name = angle.atoms[0].atomName
                a2Name = angle.atoms[1].atomName
                a3Name = angle.atoms[2].atomName
                line = "ANGLe %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name)
                topFile.write(line)

        if not autoDihFlag:
            topFile.write("\n{ Proper Dihedrals: name1 name2 name3 name4 }\n")
            for item in self.condensedProperDihedrals:
                for dih in item:
                    l = len(item)
                    a1Name = dih.atoms[0].atomName
                    a2Name = dih.atoms[1].atomName
                    a3Name = dih.atoms[2].atomName
                    a4Name = dih.atoms[3].atomName
                    line = "DIHEdral %-5s %-5s %-5s %-5s\n" % (a1Name, a2Name,
                                                               a3Name, a4Name)
                    break
                topFile.write(line)

        topFile.write("\n{ Improper Dihedrals: aName1 aName2 aName3 aName4 }\n")
        for dih in self.improperDihedrals:
            a1Name = dih.atoms[0].atomName
            a2Name = dih.atoms[1].atomName
            a3Name = dih.atoms[2].atomName
            a4Name = dih.atoms[3].atomName
            line = "IMPRoper %-5s %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name,
                                                       a4Name)
            topFile.write(line)

        topFile.write("\nEND {RESIdue %s}\n" % self.residueLabel[0])

        topFile.write("\nset echo=true end\n")

        #print "Writing CNS INP file\n"
        inpFile.write("Remarks " + head % (inp, date))
        inpData = \
"""
topology
  @%(CNS_top)s
end

parameters
  @%(CNS_par)s
  nbonds
      atom cdie shift eps=1.0  e14fac=0.4   tolerance=0.5
      cutnb=9.0 ctonnb=7.5 ctofnb=8.0
      nbxmod=5 vswitch wmin 1.0
  end
  remark dielectric constant eps set to 1.0
end

flags exclude elec ? end

segment name="    "
  chain
   coordinates @%(NEW_pdb)s
  end
end
coordinates @%(NEW_pdb)s

! Remarks If you want to shake up the coordinates a bit ...
 do (x=x+rand(10)-5) (all)
 do (y=y+rand(10)-5) (all)
 do (z=z+rand(10)-5) (all)
 write coordinates output=%(CNS_ran)s end

print threshold=0.02 bonds
print threshold=3.0 angles
print threshold=3.0 dihedrals
print threshold=3.0 impropers

Remarks Do Powell energy minimisation
minimise powell
  nstep=250 drop=40.0
end

write coordinates output=%(CNS_min)s end
write structure   output=%(CNS_psf)s end

! constraints interaction (not hydro) (not hydro) end

print threshold=0.02 bonds
print threshold=3.0 angles
print threshold=3.0 dihedrals
print threshold=3.0 impropers

flags exclude * include vdw end energy end
distance from=(not hydro) to=(not hydro) cutoff=2.6 end

stop
"""
        dictInp = {}
        dictInp['CNS_top'] = top
        dictInp['CNS_par'] = par
        dictInp['NEW_pdb'] = pdb
        dictInp['CNS_min'] = self.baseName+'_NEW_min.pdb'
        dictInp['CNS_psf'] = self.baseName+'_CNS.psf'
        dictInp['CNS_ran'] = self.baseName+'_rand.pdb'
        line = inpData % dictInp
        inpFile.write(line)

class ACTopol(AbstractTopol):
    """
        Class to build the AC topologies (Antechamber AmberTools 1.2)
    """

    def __init__(self, inputFile, chargeType = 'bcc', chargeVal = None,
            multiplicity = '1', atomType = 'gaff', force = False, basename = None,
            debug = False, outTopol = 'all', engine = 'tleap', allhdg = False,
            timeTol = 36000):

        self.debug = debug
        self.inputFile = os.path.basename(inputFile)
        self.rootDir = os.path.abspath('.')
        self.absInputFile = os.path.abspath(inputFile)
        if not os.path.exists(self.absInputFile):
            self.printWarn("input file doesn't exist")
        base, ext = os.path.splitext(self.inputFile)
        base = basename or base
        self.baseName = base # name of the input file without ext.
        self.timeTol = timeTol
        self.printDebug("Max execution time tolerance is %s" % elapsedTime(self.timeTol))
        self.ext = ext
        self.extOld = ext
        self.homeDir = self.baseName + '.acpypi'
        self.chargeType = chargeType
        self.chargeVal = chargeVal
        self.multiplicity = multiplicity
        self.atomType = atomType
        self.force = force
        self.engine = engine
        self.allhdg = allhdg
        self.acExe = getoutput('which antechamber') or '' # '/Users/alan/Programmes/antechamber-1.27/exe/antechamber'
        if not os.path.exists(self.acExe):
            self.printError("no 'antechamber' executable!")
            return None
        self.tleapExe = getoutput('which tleap') or ''
        self.sleapExe = getoutput('which sleap') or ''
        self.parmchkExe = getoutput('which parmchk') or ''
        self.babelExe = getoutput('which babel') or ''
        if not os.path.exists(self.babelExe):
            if self.ext != '.mol2':
                self.printError("no 'babel' executable... aborting!")
                self.printError("use only MOL2 file as input")
                sys.exit(1)
            else:
                self.printWarn("no 'babel' executable, no PDB file as input can be used!")
        acBase = base + '_AC'
        self.acBaseName = acBase
        self.acXyzFileName = acBase + '.inpcrd'
        self.acTopFileName = acBase + '.prmtop'
        self.acFrcmodFileName = acBase +'.frcmod'
        self.tmpDir = os.path.join(self.rootDir,'.acpypi_tmp_%s' % os.path.basename(base))
        self.setResNameCheckCoords()
        self.guessCharge()
        acMol2FileName = '%s_%s_%s.mol2' % (base, chargeType, atomType)
        self.acMol2FileName = acMol2FileName
        self.charmmBase = '%s_CHARMM' % base
        self.outTopols = [outTopol]
        if outTopol == 'all':
            self.outTopols = outTopols
        self.acParDict = {'base' : base, 'ext' : ext[1:], 'acBase': acBase,
                          'acMol2FileName' : acMol2FileName, 'res' : self.resName,
                          'leapAmberFile':leapAmberFile,
                          'leapGaffFile':leapGaffFile}

class MolTopol(ACTopol):
    """"
        Class to write topologies and parameters files for several applications

        http://amber.scripps.edu/formats.html (not updated to amber 10 yet)

        Parser, take information in AC xyz and top files and convert to objects

        INPUTS: acFileXyz and acFileTop
        RETURN: molTopol obj or None
    """
    def __init__(self, acTopolObj = None, acFileXyz = None, acFileTop = None,
                 debug = False, basename = None):

        self.allhdg = False
        self.debug = debug
        self.inputFile = acFileTop
        if acTopolObj:
            if not acFileXyz: acFileXyz = acTopolObj.acXyzFileName
            if not acFileTop: acFileTop = acTopolObj.acTopFileName
            self._parent = acTopolObj
            self.allhdg = self._parent.allhdg
            self.debug = self._parent.debug
            self.inputFile = self._parent.inputFile
        if not os.path.exists(acFileXyz) and not os.path.exists(acFileTop):
            self.printError("Files '%s' and '%s' don't exist")
            self.printError("molTopol object won't be created")
            return None

        self.xyzFileData = open(acFileXyz, 'r').readlines()
        self.topFileData = open(acFileTop, 'r').readlines()
        self.printDebug("prmtop and inpcrd files loaded")

#        self.pointers = self.getFlagData('POINTERS')

        self.getResidueLabel()
        if len(self.residueLabel) > 1:
            self.baseName = basename or os.path.splitext(os.path.basename(acFileTop))[0] # 'solute'
        else:
            self.baseName = basename or self.residueLabel[0] # 3 caps letters
        if acTopolObj:
            self.baseName = basename or acTopolObj.baseName
        self.printDebug("basename defined = '%s'" % self.baseName)

        self.getAtoms()

        self.getBonds()

        self.getAngles()

        self.getDihedrals()

        #self.setAtomPairs()

        #self.getExcludedAtoms()

        # a list of FLAGS from acTopFile that matter
#        self.flags = ( 'POINTERS', 'ATOM_NAME', 'CHARGE', 'MASS', 'ATOM_TYPE_INDEX',
#                  'NUMBER_EXCLUDED_ATOMS', 'NONBONDED_PARM_INDEX',
#                  'RESIDUE_LABEL', 'BOND_FORCE_CONSTANT', 'BOND_EQUIL_VALUE',
#                  'ANGLE_FORCE_CONSTANT', 'ANGLE_EQUIL_VALUE',
#                  'DIHEDRAL_FORCE_CONSTANT', 'DIHEDRAL_PERIODICITY',
#                  'DIHEDRAL_PHASE', 'AMBER_ATOM_TYPE' )

class Atom:
    """
        Charges in prmtop file has to be divide by 18.2223 to convert to charge
        in units of the electron charge.
        To convert ACOEF and BCOEF to r0 (Ang.) and epsilon (kcal/mol), as seen
        in gaff.dat for example; same atom type (i = j):
            r0 = 1/2 * (2 * ACOEF/BCOEF)^(1/6)
            epsilon = 1/(4 * A) * BCOEF^2
        To convert r0 and epsilon to ACOEF and BCOEF
            ACOEF = sqrt(ep_i * ep_j) * (r0_i + r0_j)^12
            BCOEF = 2 * sqrt(ep_i * ep_j) * (r0_i + r0_j)^6
                  = 2 * ACOEF/(r0_i + r0_j)^6
        where index i and j for atom types.
        Coord is given in Ang. and mass in Atomic Mass Unit.
    """
    def __init__(self, atomName, atomType, id, resid, mass, charge, coord):
        self.atomName = atomName
        self.atomType = atomType
        self.id = id
        self.resid = resid
        self.mass = mass
        self.charge = charge / qConv
        self.coords = coord

class AtomType:
    """
        AtomType per atom in gaff or amber.
    """
    def __init__(self, atomTypeName, mass, ACOEF, BCOEF):
        self.atomTypeName = atomTypeName
        self.mass = mass
        self.ACOEF = ACOEF
        self.BCOEF = BCOEF

class Bond:
    """
        attributes: pair of Atoms, spring constant (kcal/mol), dist. eq. (Ang)
    """
    def __init__(self, atoms, kBond, rEq):
        self.atoms = atoms
        self.kBond = kBond
        self.rEq = rEq

class Angle:
    """
        attributes: 3 Atoms, spring constant (kcal/mol/rad^2), angle eq. (rad)
    """
    def __init__(self, atoms, kTheta, thetaEq):
        self.atoms = atoms
        self.kTheta = kTheta
        self.thetaEq = thetaEq # rad, to convert to degree: thetaEq * 180/Pi

class Dihedral:
    """
        attributes: 4 Atoms, spring constant (kcal/mol), periodicity,
        phase (rad)
    """
    def __init__(self, atoms, kPhi, period, phase):
        self.atoms = atoms
        self.kPhi = kPhi
        self.period = period
        self.phase = phase # rad, to convert to degree: kPhi * 180/Pi

if __name__ == '__main__':
    t0 = time.time()
    print header
    argsDict = parseArgs(sys.argv[1:])
    iF = argsDict.get('-i', None)
    bn = argsDict.get('-b', None)
    to = int(argsDict.get('-s', 36000))
    cT = argsDict.get('-c', 'bcc')
    cV = argsDict.get('-n', None)
    mt = argsDict.get('-m', '1')
    at = argsDict.get('-a', 'gaff')
    ot = argsDict.get('-o', 'all')
    en = argsDict.get('-e', 'tleap')
    amb2gmx = argsDict.get('amb2gmx')
    prmtop = argsDict.get('-p', None)
    inpcrd = argsDict.get('-x', None)
    fs = False
    dg = False
    tt = False
    if '-f' in argsDict.keys(): fs = True
    if '-d' in argsDict.keys(): dg = True
    if '-t' in argsDict.keys(): tt = True

    try:
        if amb2gmx:
            print "Converting Amber input files to Gromacs ..."
            system = MolTopol(acFileXyz = inpcrd, acFileTop = prmtop, debug = dg, basename = bn)
            system.printDebug("prmtop and inpcrd files parsed")
            system.writeGromacsTopolFiles(amb2gmx = True)
        else:
            molecule = ACTopol(iF, chargeType = cT, chargeVal = cV, debug = dg,
                               multiplicity = mt, atomType = at, force = fs,
                               outTopol = ot, engine = en, allhdg=tt, basename = bn,
                               timeTol = to)

            if not molecule.acExe:
                molecule.printError("no 'antechamber' executable... aborting!")
                sys.exit(1)

            molecule.createACTopol()
            molecule.createMolTopol()
        acpypiFailed = False
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print "ACPYPI FAILED: %s" % exceptionValue
        if dg:
            traceback.print_tb(exceptionTraceback, file=sys.stdout)
        acpypiFailed = True

    execTime = int(round(time.time() - t0))
    if execTime == 0:
        msg = "less than a second"
    else:
        msg = elapsedTime(execTime)
    try: rmtree(molecule.tmpDir)
    except: pass
    print "Total time of execution: %s" % msg
    if acpypiFailed: sys.exit(1)


