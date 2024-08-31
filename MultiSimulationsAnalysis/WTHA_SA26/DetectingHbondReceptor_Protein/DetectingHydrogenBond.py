import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)
import sys
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import pool, Manager
import pandas as pd
import json
import re


def HydrogenBondDetection(Environ,frame_index,ResiID, SEGID):

	print(frame_index)
	Environ.universe.trajectory[frame_index]

	# Hbond donor and acceptor selection
	hbonds = HBA(universe=Environ,between=[f"resid {ResiID} and segid {SEGID}",f"not (resid {ResiID} and segid {SEGID})"])
	
	protein_hydrogens_sel = hbonds.guess_hydrogens(f"not (resid {ResiID} and segid {SEGID})")
	protein_acceptors_sel = hbonds.guess_acceptors(f"not (resid {ResiID} and segid {SEGID})")	

	resid_hydrogens_sel = hbonds.guess_hydrogens(f"resid {ResiID} and segid {SEGID}")
	resid_acceptors_sel = hbonds.guess_acceptors(f"resid {ResiID} and segid {SEGID}")

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
	# Target Residue chain ID for Hbond contact detection
	ChainID=sys.argv[5]


	# Load topology and trajectory to MDAnalysis
	u = MDAnalysis.Universe(TPR, Traj)

	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	uPDB = MDAnalysis.Universe(PDBFile)

	# Converting PDB to TPR
	for atomTPR, atomPDB in zip(u.atoms, uPDB.atoms):
		if (atomPDB.resid == int(ResID)) and (atomPDB.segid == ChainID):
			NewResID=atomTPR.resid
			NewChainID=atomTPR.segid
		
	# Initialize Dict to Write
	DictList=Manager().dict()

	# Multiprocessing
	pool = mp.Pool(mp.cpu_count())
	pool.starmap(HydrogenBondDetection, [(u,ts,NewResID,NewChainID) for ts in range(len(u.universe.trajectory))])

	# Write to a file
	FileName=f"HydrogenBondDetection.json"

	with open(FileName, "w") as outfile:
		json.dump(DictList.copy(), outfile)























