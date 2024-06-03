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
    # ANE5
    
    ANE5_HydrogenBondPerProteinRes={}	
    for ProteinRes in np.unique([atom.resid for atom in uPDB.select_atoms("protein")]):
        ANE5_HydrogenBondPerProteinRes.update({str(ProteinRes):0})

    #BGAL
    
    BGAL_HydrogenBondPerProteinRes={}	
    for ProteinRes in np.unique([atom.resid for atom in uPDB.select_atoms("protein")]):
        BGAL_HydrogenBondPerProteinRes.update({str(ProteinRes):0})
        
        
    #AGLC
    
    AGLC_HydrogenBondPerProteinRes={}	
    for ProteinRes in np.unique([atom.resid for atom in uPDB.select_atoms("protein")]):
        AGLC_HydrogenBondPerProteinRes.update({str(ProteinRes):0})
    

    # Read HydrogenFile	
    with open(HydrogenFile) as f:
        HydrogenFileContent = f.read().splitlines()
    
  
    # Read each hydrogen bond
    for line in HydrogenFileContent:
        HydrogenBondANE5=[]
        HydrogenBondBGAL=[]
        HydrogenBondAGLC=[]
                     
        if "ANE5" in np.unique([bond.split("_")[1] for bond in line.split()]):

           for bond in line.split():
               if "ANE5" not in bond:
                   HydrogenBondANE5.append(bond)
                   
           ProteinResHydrogenBond_ANE5=np.unique([element.split("_")[0] for element in HydrogenBondANE5])[0]
           
           # Update ANE5_HydrogenBondPerProteinRes
           ANE5_HydrogenBondPerProteinRes[ProteinResHydrogenBond_ANE5]+=1
	
        elif "BGAL" in np.unique([bond.split("_")[1] for bond in line.split()]):

           for bond in line.split():
               if "BGAL" not in bond:
                   HydrogenBondBGAL.append(bond)
                   
           ProteinResHydrogenBond_BGAL=np.unique([element.split("_")[0] for element in HydrogenBondBGAL])[0]
           
           # Update BGAL_HydrogenBondPerProteinRes
           BGAL_HydrogenBondPerProteinRes[ProteinResHydrogenBond_BGAL]+=1


        elif "AGLC" in np.unique([bond.split("_")[1] for bond in line.split()]):
                   
           for bond in line.split():
               if "AGLC" not in bond:
                   HydrogenBondAGLC.append(bond)
                   
           ProteinResHydrogenBond_AGLC=np.unique([element.split("_")[0] for element in HydrogenBondAGLC])[0]
           
           # Update AGLC_HydrogenBondPerProteinRes
           AGLC_HydrogenBondPerProteinRes[ProteinResHydrogenBond_AGLC]+=1


    # Write to JsonFile	
    with open("ANE5_HydrogenBondCountDict.json", "w") as outfile:
        json.dump(ANE5_HydrogenBondPerProteinRes, outfile)

    with open("BGAL_HydrogenBondCountDict.json", "w") as outfile:
        json.dump(BGAL_HydrogenBondPerProteinRes, outfile)

    with open("AGLC_HydrogenBondCountDict.json", "w") as outfile:
        json.dump(AGLC_HydrogenBondPerProteinRes, outfile)


