#!/usr/bin/python

import mmtools.modellertools.modelPDB as modelPDB
import sys

files = [["1YRI.pdb", "1YRI_wtv.seq"],
    ["2F4K.pdb", "2F4K_supv.seq"],
    ["1W4E.pdb", "1W4E_psbd.seq"],
    ["1BAL.pdb", "1BAL_bbl.seq"],
    ["1I6C.pdb", "1I6C_ww.seq"],
    ["1I6C.pdb", "1I6C_wwmut.seq"],
    ["1E0L.pdb", "1E0L_wwmouse.seq"],
    ["1DIV.pdb", "1DIV_ntl9.seq"],
    ["1DIV.pdb", "1DIV_ntl9mut.seq"], # check with vince about sequence/struct
    ["1ENH.pdb", "1ENH_enhd.seq"],
    ["1IDY.pdb", "1IDY_cmyb.seq"],
    ["1IDY.pdb", "1IDY_dP174.seq"],
    ["1BDD.pdb", "1BDD_proA.seq"], # check with xuhui about protein A
    ["2GB1.pdb", "2GB2_proG.seq"],
    ["2PTL.pdb", "2PTL_proL.seq"]
]

for [pdb_file, seq_file] in files:
    print "PDB File: " + pdb_file
    print "Sequence File: " + seq_file

    ext_ind = seq_file.rindex(".")
    pdbout = seq_file[0:ext_ind] + ".pdb"

    myModel = modelPDB.ModelPDB()
    myModel.makeModel(pdb_file, seq_file, pdbout)


