import json
import pandas as pd
import sys
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)

if __name__ == "__main__":


	# "HydrogenBondDetection.json" file, output of the DetectingHydrogenBond.py program
	HydrogenFile=sys.argv[1]
	# Input the topology in the form of .tpr file
	TPR=sys.argv[2]	
	# Input the topology in the form of .pdb file
	PDBFile=sys.argv[3]

	# Atom information in universe u using .tpr file as topology
	# Load topology and trajectory to MDAnalysis
	uTPR = MDAnalysis.Universe(TPR)

	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	uPDB = MDAnalysis.Universe(PDBFile)

	# Converting TPR to PDB
	
	TPR_PDB_data=[]
	for atomTPR, atomPDB in zip(uTPR.atoms, uPDB.atoms):
		TPR_PDB_data.append([int(atomTPR.id), atomTPR.name, atomTPR.resname, atomTPR.segid, atomTPR.resid,
				int(atomPDB.id), atomPDB.name, atomPDB.resname, atomPDB.segid, atomPDB.resid])

	TPR_PDB_df=pd.DataFrame(TPR_PDB_data, columns=["AtomTPR","AtomTypeTPR","ResNameTPR","ChainIDTPR","ResIDTPR",
						"AtomPDB","AtomTypePDB","ResNamePDB","ChainIDPDB","ResIDPDB"])  
			
	# Opening Hydrogen bondJSON file
	f=open(HydrogenFile)
	data=json.load(f)
	
	# Read the hydrogen bond dict
	# Convert TPR to PDB
	for key, value in data.items():
		
		for bond in value:
			TPRDonorID=int(bond[1]) 
			TPRHydrogenID=int(bond[2]) 
			TPRAcceptorID=int(bond[3])
			
			DonorIDidx=TPR_PDB_df.index[TPR_PDB_df["AtomTPR"]==TPRDonorID].tolist()[0]
			HydrogenIDidx=TPR_PDB_df.index[TPR_PDB_df["AtomTPR"]==TPRHydrogenID].tolist()[0]
			AcceptorIDidx=TPR_PDB_df.index[TPR_PDB_df["AtomTPR"]==TPRAcceptorID].tolist()[0]

			# Donor in pdb format
			PDBDonorIDRes1=TPR_PDB_df["ResIDPDB"][DonorIDidx]
			PDBDonorIDRes2=TPR_PDB_df["ChainIDPDB"][DonorIDidx]
			PDBDonorIDRes3=TPR_PDB_df["ResNamePDB"][DonorIDidx]
			PDBDonorIDRes4=TPR_PDB_df["AtomTypePDB"][DonorIDidx]
			DonorIDRes=f"{PDBDonorIDRes1}{PDBDonorIDRes2}_{PDBDonorIDRes3}_{PDBDonorIDRes4}"
			
			# Hydrogen
			PDBHydrogenIDRes1=TPR_PDB_df["ResIDPDB"][HydrogenIDidx]
			PDBHydrogenIDRes2=TPR_PDB_df["ChainIDPDB"][HydrogenIDidx]
			PDBHydrogenIDRes3=TPR_PDB_df["ResNamePDB"][HydrogenIDidx]
			PDBHydrogenIDRes4=TPR_PDB_df["AtomTypePDB"][HydrogenIDidx]
			HydrogenIDRes=f"{PDBHydrogenIDRes1}{PDBHydrogenIDRes2}_{PDBHydrogenIDRes3}_{PDBHydrogenIDRes4}"			


			# Acceptor
			PDBAcceptorIDRes1=TPR_PDB_df["ResIDPDB"][AcceptorIDidx]
			PDBAcceptorIDRes2=TPR_PDB_df["ChainIDPDB"][AcceptorIDidx]
			PDBAcceptorIDRes3=TPR_PDB_df["ResNamePDB"][AcceptorIDidx]
			PDBAcceptorIDRes4=TPR_PDB_df["AtomTypePDB"][AcceptorIDidx]
			AcceptorIDRes=f"{PDBAcceptorIDRes1}{PDBAcceptorIDRes2}_{PDBAcceptorIDRes3}_{PDBAcceptorIDRes4}"
			
			
			# Write to a file
			OutFile=open("HydrogenBondPerTimeStep.txt", "a")
			OutFile.write(f"{DonorIDRes} {HydrogenIDRes} {AcceptorIDRes}\n")
			OutFile.close()
       			
			


            
            
            
            
