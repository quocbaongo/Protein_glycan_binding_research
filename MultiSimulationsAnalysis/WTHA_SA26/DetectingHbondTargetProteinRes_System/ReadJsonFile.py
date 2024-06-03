import json
import pandas as pd
import sys
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)

if __name__ == "__main__":

    # "HydrogenBondDetection.json" file, output of the DetectingHydrogenBond.py program
    HydrogenFile=sys.argv[1]
    # Input trajectory in the form of .xtc file
    Traj=sys.argv[2]
    # Input the topology in the form of .tpr file
    TPR=sys.argv[3]
    # Input the topology in the form of .pdb file
    PDBFile=sys.argv[4]
    
    # Load topology and trajectory to MDAnalysis
    u = MDAnalysis.Universe(TPR, Traj)
    
    # Residue ID in universe u is incorrect compared to the input struct
    # Using the input PDBFile to map them correctly

    # Atom information in universe u using .tpr file as topology
    TPRData=[]
    for atom in u.select_atoms("all"):
    	#print(atom)
    	TPRData.append([int(atom.id),atom.name,atom.resname,atom.segid,atom.resid])
    
    TPR_df=pd.DataFrame(TPRData, columns=["Atom","AtomType","ResName","ChainID","ResID"])  	
    

    # Atom information in universe u using .pdb file as topology
    # Load topology and trajectory to MDAnalysis
    uPDB = MDAnalysis.Universe(PDBFile)
    
    PDBData=[]
    for atom in uPDB.select_atoms("all"):

    	PDBData.append([int(atom.id),atom.name,atom.resname,atom.segid,atom.resid])
    
    PDB_df=pd.DataFrame(PDBData, columns=["Atom","AtomType","ResName","ChainID","ResID"])
    
    # TPR to PDB converter
    TPRToPDB={}
    
    for i,j in zip(PDB_df.index,TPR_df.index):
        # AtomID in TPR: [AtomID in PDB, ResID in PDB that corresponds to the on in TPR, ResName in PDB that corresponds to the on in TPR, Atom type in PDB]
        TPRToPDB.update({TPR_df["Atom"][j]:[PDB_df["Atom"][i], PDB_df["ResID"][i], PDB_df["ResName"][i], PDB_df["AtomType"][i]]})

    #print(TPRToPDB)

    # Opening Hydrogen bondJSON file
    f=open(HydrogenFile)
    data=json.load(f)

    # Read the hydrogen bond dict
    CountProgress=0
    
    for key, value in data.items():
        print(CountProgress)
        CountProgress+=1
        for bond in value:
            
            DonorID=int(bond[1])
            DonorIDRes=f"{TPRToPDB[DonorID][1]}_{TPRToPDB[DonorID][2]}_{TPRToPDB[DonorID][3]}"
           
            HydrogenID=int(bond[2])
            HydrogenIDRes=f"{TPRToPDB[HydrogenID][1]}_{TPRToPDB[HydrogenID][2]}_{TPRToPDB[HydrogenID][3]}"
            
            AcceptorID=int(bond[3])
            AcceptorIDRes=f"{TPRToPDB[AcceptorID][1]}_{TPRToPDB[AcceptorID][2]}_{TPRToPDB[AcceptorID][3]}"
            
            
            # Write to a file
            OutFile=open("HydrogenBondPerTimeStep.txt", "a")
            OutFile.write(f"{DonorIDRes} {HydrogenIDRes} {AcceptorIDRes}\n")	
            OutFile.close()

