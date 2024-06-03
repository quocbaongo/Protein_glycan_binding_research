import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (HydrogenBondAnalysis as HBA)
import sys
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import pool, Manager
import json
def HydrogenBondDetection(Environ,frame_index):

	print(frame_index)	
	Environ.universe.trajectory[frame_index]

	# Hbond donor and acceptor selection
	hbonds = HBA(universe=Environ,between=["resname AGLCNA or resname BGAL or resname ANE5AC", "protein"])
	protein_hydrogens_sel = hbonds.guess_hydrogens("protein")
	protein_acceptors_sel = hbonds.guess_acceptors("protein")

	LIG_hydrogens_sel = hbonds.guess_hydrogens("resname AGLCNA or resname BGAL or resname ANE5AC")
	LIG_acceptors_sel = hbonds.guess_acceptors("resname AGLCNA or resname BGAL or resname ANE5AC")	

	hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({LIG_hydrogens_sel})"
	hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({LIG_acceptors_sel})"		
    
	hbonds.run(start=int(frame_index),stop=(int(frame_index)+1))

	# Write to list
	DictList.update({str(frame_index):hbonds.results.hbonds.tolist()})

if __name__ == "__main__":

	# Input trajectory in the form of .xtc file
	TRAJ=sys.argv[1]
	
	# Input the topology in the form of .tpr file
	TPR=sys.argv[2]

	# Initiate MDAnalysis
	u=MDAnalysis.Universe(TPR, TRAJ)

	# Initialize Dict to Write
	DictList=Manager().dict()	

	# Multiprocessing
	pool = mp.Pool(mp.cpu_count()) 
	pool.starmap(HydrogenBondDetection, [(u,ts) for ts in range(len(u.universe.trajectory))]) 	

	# Write to a file
	FileName="HydrogenBondDetection.json"
			
	with open(FileName, "w") as outfile:
		json.dump(DictList.copy(), outfile)

