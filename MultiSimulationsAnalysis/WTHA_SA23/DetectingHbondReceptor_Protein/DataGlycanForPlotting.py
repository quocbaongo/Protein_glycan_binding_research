import json
import pandas as pd
import sys
import MDAnalysis
import numpy as np
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)

AminoAcid=["VAL","ILE","LEU","GLU","GLN","ASP","ASN","HIS","TRP","PHE","TYR","ARG","LYS",
        "SER","THR","MET","ALA","GLY","PRO","CYS"]

if __name__ == "__main__":

    # "HydrogenBondPerTimeStep.txt" file, output of the ReadJsonFile.py program
    HydrogenFile=sys.argv[1]
    # Input trajectory in the form of .xtc file
    Traj=sys.argv[2]
    # Input the topology in the form of .pdb file
    PDBFile=sys.argv[3]

    # Atom information in universe u using .pdb file as topology
    # Load topology and trajectory to MDAnalysis
    uPDB = MDAnalysis.Universe(PDBFile)
    
    # Create dictionary for counting hydrogen bond
    HydrogenBondPerProteinRes={}	
    for ProteinRes in np.unique([atom.resid for atom in uPDB.select_atoms("protein")]):
        HydrogenBondPerProteinRes.update({str(ProteinRes):0})

    # Read HydrogenFile	
    with open(HydrogenFile) as f:
        HydrogenFileContent = f.read().splitlines()	    
        
     # Read each hydrogen bond
    for line in HydrogenFileContent:
        HydrogenBondConnection=[]
        for atom in line.split():
            if atom.split("_")[1] in AminoAcid:
                HydrogenBondConnection.append(atom.split("_")[0])
            
        ProteinResHydrogenBond=np.unique(HydrogenBondConnection)[0]        

        # Update HydrogenBondPerProteinRes
        HydrogenBondPerProteinRes[ProteinResHydrogenBond]+=1

    # Write to JsonFile			
    with open("HydrogenBond_WholeGlycan_CountDict.json", "w") as outfile:
        json.dump(HydrogenBondPerProteinRes, outfile)

