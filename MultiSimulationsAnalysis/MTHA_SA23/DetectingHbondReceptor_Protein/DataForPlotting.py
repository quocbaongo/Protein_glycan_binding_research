import json
import pandas as pd
import sys
import MDAnalysis
import numpy as np
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)

AminoAcid=["VAL","ILE","LEU","GLU","GLN","ASP","ASN","HIS",
	"TRP","PHE","TYR","ARG","LYS","SER","THR","MET",
	"ALA","GLY","PRO","CYS"]

if __name__ == "__main__":
	
	# "HydrogenBondPerTimeStep.txt" file, output of the ReadJsonFile.py program
	HydrogenFile=sys.argv[1]	
	# Input the topology in the form of .pdb file
	PDBFile=sys.argv[2]
	# Target Residue ID for Hbond contact detection
	ResID=sys.argv[3]	
	# Chain ID for Hbond contact detection
	ChainID=sys.argv[4]
	
	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	uPDB = MDAnalysis.Universe(PDBFile)
    
	# Create dictionary for counting hydrogen bond
	HydrogenBondPerProteinRes={}
	for AA in uPDB.select_atoms("all"):
		HydrogenBondPerProteinRes.update({f"{AA.resid}{AA.segid}": 0})

    
	# Read HydrogenFile
	with open(HydrogenFile) as f:
		HydrogenFileContent = f.read().splitlines()    
    
    
	# Read each hydrogen bond
	for line in HydrogenFileContent:
		HydrogenBondConnection=[]
		
		for atom in line.split():
			HydrogenBondConnection.append(atom.split("_")[0])
		
		HydrogenBondConnection=list(dict.fromkeys(HydrogenBondConnection))
		
		for element in HydrogenBondConnection:
			if element != str(f"{ResID}{ChainID}"):
				HydrogenBondPerProteinRes[element] += 1
	
	# Write to JsonFile
	with open(f"HydrogenBondCountDict.json", "w") as outfile:
		json.dump(HydrogenBondPerProteinRes, outfile)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    	
