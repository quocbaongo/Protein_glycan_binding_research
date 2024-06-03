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
    # Target Residue ID for Hbond contact detection
    ResID=sys.argv[4]
    
    # Atom information in universe u using .pdb file as topology
    # Load topology and trajectory to MDAnalysis
    uPDB = MDAnalysis.Universe(PDBFile)
    
    # Create dictionary for counting hydrogen bond
    HydrogenBondPerProteinRes={}	
    for ProteinRes in np.unique([atom.resid for atom in uPDB.select_atoms("protein")]):
        HydrogenBondPerProteinRes.update({str(ProteinRes):0})

    HydrogenBondPerProteinRes.update({"Glycan":0})

    # Read HydrogenFile	
    with open(HydrogenFile) as f:
        HydrogenFileContent = f.read().splitlines()	    
        
     # Read each hydrogen bond
    for line in HydrogenFileContent:
        HydrogenProteinBondConnection=[]
        HydrogenGlycanBondConnection=[]
        
        for atom in line.split():
            if int(atom.split("_")[0]) != int(ResID) and atom.split("_")[1] in AminoAcid:
                HydrogenProteinBondConnection.append(atom.split("_")[0])

            elif int(atom.split("_")[0]) != int(ResID) and atom.split("_")[1] in ["ANE5","BGLC","BGAL"]:
                HydrogenGlycanBondConnection.append(atom.split("_")[0])
                
        if len(HydrogenProteinBondConnection) != 0:
        	ProteinResHydrogenBond=np.unique(HydrogenProteinBondConnection)[0]        

		# Update HydrogenBondPerProteinRes
        	HydrogenBondPerProteinRes[ProteinResHydrogenBond]+=1

        if len(HydrogenGlycanBondConnection) != 0:
        	ProteinResHydrogenBond=np.unique(HydrogenGlycanBondConnection)[0]
        	
        	# Update HydrogenBondPerProteinRes
        	HydrogenBondPerProteinRes["Glycan"]+=1
        
    # Write to JsonFile			
    with open(f"HydrogenBondCountDict{ResID}.json", "w") as outfile:
        json.dump(HydrogenBondPerProteinRes, outfile)

