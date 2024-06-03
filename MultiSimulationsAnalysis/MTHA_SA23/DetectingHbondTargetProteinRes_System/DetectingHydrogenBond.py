import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)
import sys
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import pool, Manager
import pandas as pd
import json

def HydrogenBondDetection(Environ,frame_index,ResiID):

	print(frame_index)	
	Environ.universe.trajectory[frame_index]

	# Hbond donor and acceptor selection
	hbonds = HBA(universe=Environ,between=[f"resid {ResiID} and protein",f"(protein and not (resid {ResiID} and protein)) or resname BGLCNA or resname BGAL or resname ANE5AC"])
	
	protein_hydrogens_sel = hbonds.guess_hydrogens(f"(protein and not (resid {ResiID} and protein)) or resname BGLCNA or resname BGAL or resname ANE5AC")
	protein_acceptors_sel = hbonds.guess_acceptors(f"(protein and not (resid {ResiID} and protein)) or resname BGLCNA or resname BGAL or resname ANE5AC")

	resid_hydrogens_sel = hbonds.guess_hydrogens(f"resid {ResiID} and protein")
	resid_acceptors_sel = hbonds.guess_acceptors(f"resid {ResiID} and protein")
	
	hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({resid_hydrogens_sel})"
	hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({resid_acceptors_sel})"
	
	hbonds.run(start=int(frame_index),stop=(int(frame_index)+1))

	# Write to list
	DictList.update({str(frame_index):hbonds.results.hbonds.tolist()})



if __name__ == "__main__":

	# Input trajectory in the form of .xtc file
	Traj=sys.argv[1]
	# Input the topology in the form of .tpr file
	TPR=sys.argv[2]
	# Input the topology in the form of .pdb file
	PDBFile=sys.argv[3]
	# Target Residue ID for Hbond contact detection
	ResID=sys.argv[4]
	
	
	# Load topology and trajectory to MDAnalysis
	u = MDAnalysis.Universe(TPR, Traj)
	
	# Residue ID in universe u is incorrect compared to the input struct
	# Using the input PDBFile to map them correctly
	
	# Atom information in universe u using .tpr file as topology
	TPRData=[]
	
	for atom in u.select_atoms("protein"):
		#print(atom)
		TPRData.append([int(atom.id),atom.type,atom.resname,atom.segid,atom.resid])
		
	TPR_df=pd.DataFrame(TPRData, columns=["Atom","AtomType","ResName","ChainID","ResID"])
	
	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	
	uPDB = MDAnalysis.Universe(PDBFile)
	
	PDBData=[]
	for atom in uPDB.select_atoms("protein"):
		PDBData.append([int(atom.id),atom.type,atom.resname,atom.segid,atom.resid])
		
	PDB_df=pd.DataFrame(PDBData, columns=["Atom","AtomType","ResName","ChainID","ResID"])
	
	# TPR to PDB converter
	TPRToPDB={}
	
	for i,j in zip(PDB_df.index,TPR_df.index):
		TPRToPDB.update({PDB_df["ResID"][i]:TPR_df["ResID"][j]})
		
	# Initialize Dict to Write
	DictList=Manager().dict()
	
	# Multiprocessing
	ResIDNew=TPRToPDB[int(ResID)]
	
	pool = mp.Pool(mp.cpu_count())
	pool.starmap(HydrogenBondDetection, [(u,ts,ResIDNew) for ts in range(len(u.universe.trajectory))])
	
	# Write to a file
	FileName=f"HydrogenBondDetection{ResID}.json"
	
	with open(FileName, "w") as outfile:
		json.dump(DictList.copy(), outfile)

