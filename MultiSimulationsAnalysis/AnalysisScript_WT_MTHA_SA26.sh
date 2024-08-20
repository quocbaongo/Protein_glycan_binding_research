GMX=gmx_2021.2

# Input directory
WorkDir=/path/to/Working/Directory

# The very first step of the analysis is to fix the periodic boundary condition of the generated trajectory
FixPBCTrajDir=/path/to/directory/containing/all/the/simulation/trajectory

# Protein-Glycan initial topology
EMTPR=/path/to/system/topology/in/tpr

# Starting conformations
StartingPDBDir=/path/to/10/Initial/ProteinGlycan/structure

set -e


# We generated an ensemble of conformations between the protein and glycan
#by concatenating all the 10 individual trajectories to 1 trajectory only.
#All the following analysis would be conducted on the concatenated trajectory.

cd $WorkDir/ConcatenatedTraj

$GMX trjcat -f $FixPBCTrajDir/frame_*/FixPBC/md_whole_cluster_nojump_fit.xtc -settime -cat <<EOF
l
c
c
c
c
c
c
c
c
c
EOF

echo "Protein_LIG" | $GMX trjconv -f $WorkDir/ConcatenatedTraj/trajout.xtc -s $WorkDir/ConcatenatedTraj/tpr.tpr -dump 0 -o $WorkDir/ConcatenatedTraj/frame0.pdb -n $WorkDir/ConcatenatedTraj/index.ndx



	#################################################################################################################################################

# The objective of this analysis is to process the concatenated trajectory "trajout.xtc" 
# to exclude all the frames, in which the glycan diffused away its binding pocket in the HA protein 

cd $WorkDir/ConcatenatedTraj

echo 'import MDAnalysis
import sys

def MeasuringDistance(Ref, Target):

	x1=Ref[0]
	y1=Ref[1]
	z1=Ref[2]
				 
	x2=Target[0]
	y2=Target[1]
	z2=Target[2]
	
	dist=((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
		
	return dist

if __name__ == "__main__":

	# Input trajectory in the form of .xtc file
	TRAJ=sys.argv[1]
	
	# Input the topology in the form of .pdb file
	TPR=sys.argv[2]

	# Distance dictionary
	DictList={}

	# Initiate MDAnalysis
	u = MDAnalysis.Universe(TPR, TRAJ)
	
	# List for Unbound frame
	UnboundFrameIdx=[]
	# Extracting protein glycan bound state
	count=0
	with MDAnalysis.Writer("BoundTrajout.xtc", multiframe=True, n_atoms=u.atoms.n_atoms) as W:
	    for ts in u.trajectory:
	        print(count)
	      
	        # Residue that belong to binding pocket
	        BindingPocket=u.select_atoms("(resid 91 or resid 92) and segid A")
	        BindingPocket+=u.select_atoms("resid 126:136 and segid A")
	        BindingPocket+=u.select_atoms("resid 141:143 and segid A")
	        BindingPocket+=u.select_atoms("resid 149:152 and segid A")
	        BindingPocket+=u.select_atoms("resid 179:191 and segid A")
	        BindingPocket+=u.select_atoms("resid 216:226 and segid A")
	        
	        BindingPocketResID=[atom.resid for atom in BindingPocket]
	        BindingPocketResID=list(set(BindingPocketResID))
	        
	        # Binding pocket center of geometry
	        BindingPocketCOM=BindingPocket.center_of_geometry()
	        
	        # Protein residue that is 1.2 nm around the glycan	
	        AroundGlycan=[atom.resid for atom in (u.select_atoms("around 12.0 (resname AGLC or resname BGAL or resname ANE5)",periodic=False))]
	        AroundGlycan=list(set(AroundGlycan))
	        print(len(AroundGlycan))
	        
	        # Receptor center of geometry
	        ReceptorHeavyAtomsCoord=u.select_atoms("(resname AGLC or resname BGAL or resname ANE5) and not name H*").center_of_geometry()
	        
	        distance=MeasuringDistance(ReceptorHeavyAtomsCoord, BindingPocketCOM)
	        
	        check=any(item in AroundGlycan for item in BindingPocketResID)
	        
	        if check == True and distance <= 10:
	            count +=1
	            W.write(u.atoms)
	        elif len(AroundGlycan) == 0:
	            UnboundFrameIdx.append(ts.time)
	        
	       	        
	# Extracting protein glycan unbound state
	with MDAnalysis.Writer("UnBoundTrajout.xtc", multiframe=True, n_atoms=u.atoms.n_atoms) as U:
	    for ts in u.trajectory:
	        if ts.time in UnboundFrameIdx:
	            U.write(u.atoms)
' > $WorkDir/ConcatenatedTraj/PostProcessingTrajout.py

python3 $WorkDir/ConcatenatedTraj/PostProcessingTrajout.py $WorkDir/ConcatenatedTraj/trajout.xtc $WorkDir/ConcatenatedTraj/frame0.pdb



	#######################################################################################



# We measured the distance between the center of geometry of the glycan and 
#center of geometry of a group of CA atoms of protein residues that 
#together form glycan binding pocket over the course of every individual simulation

mkdir $WorkDir/ProteinReceptorDistance
cd $WorkDir/ProteinReceptorDistance

# Generate a python file to measure the distance between ligand and protein
#over the course of the simulation. Through visual inspection using PyMol, we
#noticed that the receptor binding pocket on the surface of the HA protein is 
#formed by a group of residue including resid 91, 92, 126-136, 141-143, 149-152
#179-191, 216-226. 
# The output of the python program is a json that indicates the distance between 
#center of geometry of the protein ligand binding pocket and center of geometry
#of ligand.

echo 'import MDAnalysis
import sys
import json

def MeasuringDistance(Ref, Target):

	x1=Ref[0]
	y1=Ref[1]
	z1=Ref[2]
				 
	x2=Target[0]
	y2=Target[1]
	z2=Target[2]
	
	dist=((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
		
	return dist

if __name__ == "__main__":

	# Input trajectory in the form of .xtc file
	TRAJ=sys.argv[1]
	
	# Input the topology in the form of .pdb file
	TPR=sys.argv[2]

	# Distance dictionary
	DictList={}

	# Initiate MDAnalysis
	u = MDAnalysis.Universe(TPR, TRAJ)

	for ts in range(len(u.universe.trajectory)):
		print(ts)	
		u.universe.trajectory[ts]
			
		BindingPocket=u.select_atoms("(resid 91 or resid 92) and segid A")		
		BindingPocket+=u.select_atoms("resid 126:136 and segid A")
		BindingPocket+=u.select_atoms("resid 141:143 and segid A")
		BindingPocket+=u.select_atoms("resid 149:152 and segid A")
		BindingPocket+=u.select_atoms("resid 179:191 and segid A")
		BindingPocket+=u.select_atoms("resid 216:226 and segid A")

		BindingPocketCOM=BindingPocket.center_of_geometry()

		# Receptor center of geometry
		ReceptorHeavyAtomsCoord=u.select_atoms("(resname AGLC or resname BGAL or resname ANE5) and not name H*").center_of_geometry()
		
		DictList.update({str(ts):MeasuringDistance(ReceptorHeavyAtomsCoord, BindingPocketCOM)})

	# Write to a file
	FileName="ReceptorCOMToBindingPocketCOM.json"
			
	with open(FileName, "w") as outfile:
		json.dump(DictList, outfile)	
' > $WorkDir/ProteinReceptorDistance/Protein_Receptor_Distance.py

for i in $FixPBCTrajDir/frame_*
do
	f="$(basename -- $i)"
	mkdir $WorkDir/ProteinReceptorDistance/$f
	cd $WorkDir/ProteinReceptorDistance/$f

	python3 $WorkDir/ProteinReceptorDistance/Protein_Receptor_Distance.py $i/FixPBC/md_whole_cluster_nojump_fit.xtc $i/FixPBC/frame0.pdb
		
done



	#######################################################################################




# Detecting the Hbond contact involving a specific residue in the protein structure
#in each frame of the concatenated trajectory
mkdir $WorkDir/DetectingHbondTargetProteinRes_System
cd $WorkDir/DetectingHbondTargetProteinRes_System

# Generate a python file to detect the Hbond contact involving a specific residue in the protein structure
# The file allows multi-cpu usage
# The output of the python file is a json file named 'HydrogenBondDetection$ResIDOfInterested.json'
# Note: in the protein-receptor complex structure, the residue name of the receptor is 'LIG'

echo 'import MDAnalysis
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
	hbonds = HBA(universe=Environ,between=[f"resid {ResiID} and protein",f"(protein and not (resid {ResiID} and protein)) or resname AGLCNA or resname BGAL or resname ANE5AC"])
	
	protein_hydrogens_sel = hbonds.guess_hydrogens(f"(protein and not (resid {ResiID} and protein)) or resname AGLCNA or resname BGAL or resname ANE5AC")
	protein_acceptors_sel = hbonds.guess_acceptors(f"(protein and not (resid {ResiID} and protein)) or resname AGLCNA or resname BGAL or resname ANE5AC")

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
		TPRData.append([int(atom.id),atom.name,atom.resname,atom.segid,atom.resid])
		
	TPR_df=pd.DataFrame(TPRData, columns=["Atom","AtomType","ResName","ChainID","ResID"])
	
	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	
	uPDB = MDAnalysis.Universe(PDBFile)
	
	PDBData=[]
	for atom in uPDB.select_atoms("protein"):
		PDBData.append([int(atom.id),atom.name,atom.resname,atom.segid,atom.resid])
		
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
' > $WorkDir/DetectingHbondTargetProteinRes_System/DetectingHydrogenBond.py


# Read and interpret the generated json file 'HydrogenBondDetection.json' by using a python file named 'ReadJsonFile.py' 
echo 'import json
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
' > $WorkDir/DetectingHbondTargetProteinRes_System/ReadJsonFile.py


# Finally we used the file "HydrogenBondPerTimeStep.txt" generated in the previous step as the input
#for the python program 'DataForPlotting.py' to create a dictionary that count the number of Hbond between the receptor
#(in the form of a small molecule) and each individual residues of the protein in the ensemble of structure
#(concatenated trajectory)
# The output of the program 'DataForPlotting.py' is a json file 'HydrogenBondCountDict.json' that is ready for plotting.

echo 'import json
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

            elif int(atom.split("_")[0]) != int(ResID) and atom.split("_")[1] in ["ANE5","AGLC","BGAL"]:
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
' > $WorkDir/DetectingHbondTargetProteinRes_System/DataForPlotting.py



for i in 91 94
do
	mkdir $WorkDir/DetectingHbondTargetProteinRes_System/Res$i
	cd $WorkDir/DetectingHbondTargetProteinRes_System/Res$i
	
	python3 $WorkDir/DetectingHbondTargetProteinRes_System/DetectingHydrogenBond.py $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/tpr.tpr $WorkDir/ConcatenatedTraj/frame0.pdb $i
	
	python3 $WorkDir/DetectingHbondTargetProteinRes_System/ReadJsonFile.py $WorkDir/DetectingHbondTargetProteinRes_System/Res$i/HydrogenBondDetection$i.json $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/tpr.tpr $WorkDir/ConcatenatedTraj/frame0.pdb
	
	python3 $WorkDir/DetectingHbondTargetProteinRes_System/DataForPlotting.py $WorkDir/DetectingHbondTargetProteinRes_System/Res$i/HydrogenBondPerTimeStep.txt $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/frame0.pdb $i

done



	#######################################################################################


# Detecting the Hbond formed between receptor (glycan) and the protein in each frame of the concatenated trajectory
mkdir $WorkDir/DetectingHbondReceptor_Protein
cd $WorkDir/DetectingHbondReceptor_Protein

# Generate a python file to detect the Hbond between the receptor and protein
# The file allows multi-cpu usage
# The output of the python file is a json file named 'HydrogenBondDetection.json'
# Note: in the protein-receptor complex structure, the residue name of the receptor is 'AGLCNA', 'BGAL' and 'ANE5AC'

echo 'import MDAnalysis
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
' > $WorkDir/DetectingHbondReceptor_Protein/DetectingHydrogenBond.py 		

# It is important to note that 
#here we use the concatenated trajectory that ONLY contain protein and ligand atoms
#Thus, we created a compatible topology that also contains ONLY protein and ligand atoms 
#(no water or ions atoms)

python3 $WorkDir/DetectingHbondReceptor_Protein/DetectingHydrogenBond.py $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/tpr.tpr


# Read and interpret the generated json file 'HydrogenBondDetection.json' by using a python file named 'ReadJsonFile.py' 
echo 'import json
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
' > $WorkDir/DetectingHbondReceptor_Protein/ReadJsonFile.py
 

python3 $WorkDir/DetectingHbondReceptor_Protein/ReadJsonFile.py $WorkDir/DetectingHbondReceptor_Protein/HydrogenBondDetection.json $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/tpr.tpr $WorkDir/ConcatenatedTraj/frame0.pdb

# Finally we used the file "HydrogenBondPerTimeStep.txt" generated in the previous step as the input
#for the python program 'DataForPlotting.py' to create a dictionary that count the number of Hbond between the receptor
#(in the form of a small molecule) and each individual residues of the protein in the ensemble of structure
#(concatenated trajectory)
# The output of the program 'DataForPlotting.py' is a json file 'HydrogenBondCountDict.json' that is ready for plotting.

echo 'import json
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
' > $WorkDir/DetectingHbondReceptor_Protein/DataGlycanForPlotting.py

python3 $WorkDir/DetectingHbondReceptor_Protein/DataGlycanForPlotting.py $WorkDir/DetectingHbondReceptor_Protein/HydrogenBondPerTimeStep.txt $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/frame0.pdb


echo 'import json
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

' > $WorkDir/DetectingHbondReceptor_Protein/DataGlycanResForPlotting.py

python3 $WorkDir/DetectingHbondReceptor_Protein/DataGlycanResForPlotting.py $WorkDir/DetectingHbondReceptor_Protein/HydrogenBondPerTimeStep.txt $WorkDir/ConcatenatedTraj/BoundTrajout.xtc $WorkDir/ConcatenatedTraj/frame0.pdb


	#######################################################################################
	
# Measuring root mean square fluctuation on all C-alpha atoms
mkdir $WorkDir/RMSF
cd $WorkDir/RMSF

echo "C-alpha" | $GMX rmsf -s $WorkDir/ConcatenatedTraj/tpr.tpr -f $WorkDir/ConcatenatedTraj/BoundTrajout.xtc -o $WorkDir/RMSF/RMSFCA.xvg	


	#######################################################################################

# Sampling the angle θ of the glycan within the HA protein binding pocket
# The glycan angle θ is formed by atom C2 of Sialic Acid (SIA) sugar, C1 of Galactose (GAL) sugar 
# and N-acetylglucosamine (GlcNac) sugar


mkdir $WorkDir/Angle
cd $WorkDir/Angle

# For complex of WT HA bound to receptor SA23, use the following commands
echo -e 'a C2 & "ANE5AC" \n a C1 & "BGAL" \n a C1 & "AGLCNA" \n  "C2_&_ANE5AC" | "C1_&_BGAL" | "C1_&_AGLCNA" \n q' | $GMX make_ndx -f $WorkDir/ConcatenatedTraj/tpr.tpr -o $WorkDir/Angle/index.ndx

echo "C2_&_ANE5AC_C1_&_BGAL_C1_&_AGLCNA" | $GMX angle -f $WorkDir/ConcatenatedTraj/BoundTrajout.xtc -n $WorkDir/Angle/index.ndx -ov $WorkDir/Angle/angleBound.xvg

echo "C2_&_ANE5AC_C1_&_BGAL_C1_&_AGLCNA" | $GMX angle -f $WorkDir/ConcatenatedTraj/trajout.xtc -n $WorkDir/Angle/index.ndx -ov $WorkDir/Angle/angleAll.xvg

# Measuring angle in each starting strct file

for((i=0;i<10;i++))
do
	echo "C2_&_ANE5AC_C1_&_BGAL_C1_&_AGLCNA" | $GMX angle -f $StartingPDBDir/frame_$i.pdb -n $WorkDir/Angle/index.ndx -ov $WorkDir/Angle/angle_start$i.xvg
done	


	#################################################################################################################################################


# Sampling the RSMD of the glycan within the HA protein binding pocket

mkdir $WorkDir/ReceptorRMSD
cd $WorkDir/ReceptorRMSD

# Create a group of receptor heavy atom
echo -e '"AGLCNA" | "BGAL" | "ANE5AC" & !a H* \n q' | $GMX make_ndx -f $WorkDir/ConcatenatedTraj/tpr.tpr -o $WorkDir/ReceptorRMSD/index.ndx

# Sampling the RMSD of the receptor heavy atoms
echo -e 'AGLCNA_BGAL_ANE5AC_&_!H*' 'AGLCNA_BGAL_ANE5AC_&_!H*' | $GMX rms -s $EMTPR -f $WorkDir/ConcatenatedTraj/BoundTrajout.xtc -n $WorkDir/ReceptorRMSD/index.ndx -o $WorkDir/ReceptorRMSD/rmsdBound.xvg

echo -e 'AGLCNA_BGAL_ANE5AC_&_!H*' 'AGLCNA_BGAL_ANE5AC_&_!H*' | $GMX rms -s $EMTPR -f $WorkDir/ConcatenatedTraj/trajout.xtc -n $WorkDir/ReceptorRMSD/index.ndx -o $WorkDir/ReceptorRMSD/rmsdAll.xvg

# Measuring starting strct RMSD
for((i=0;i<10;i++))
do
	echo -e 'AGLCNA_BGAL_ANE5AC_&_!H*' 'AGLCNA_BGAL_ANE5AC_&_!H*' | $GMX rms -s $EMTPR -f $StartingPDBDir/frame_$i.pdb -n $WorkDir/ReceptorRMSD/index.ndx -o $WorkDir/ReceptorRMSD/rmsd_start$i.xvg
done


	#######################################################################################################################################################

mkdir $WorkDir/RBS_volume_calculation
cd $WorkDir/RBS_volume_calculation

# Generate a python file to measure the protein receptor binding pocket
#over the course of the simulation

echo 'import numpy as np
from scipy.spatial import ConvexHull
import MDAnalysis
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import Pool
from functools import partial
import sys

def MeasuringVolume(Environ, frame_index, SelectedRes, FileOut):

	
	Environ.universe.trajectory[frame_index]

	CoordGeometry = []
			
	for element in SelectedRes:
		for i in range(element[0], element[1]+1):
			CoordGeometry.append(Environ.select_atoms(f'resid {i} and segid {element[2]}').positions[0])
			
	CoordGeometry = np.array(CoordGeometry)
	hull = ConvexHull(CoordGeometry)
			
	# Write in file
	Volume = open(FileOut, 'a')
	Volume.write(f'{frame_index} {hull.volume}\n')
	Volume.close()


if __name__ == '__main__':

	# Input trajectory in the form of .xtc file
	TRAJ=sys.argv[1]
	
	# Input the topology in the form of .pdb file
	TPR=sys.argv[2]
	
	# File in .txt format to define amino acids that belongs to the target region
	ResIDRegion=sys.argv[3]
	
	# Initiate universe
	u = MDAnalysis.Universe(TPR, TRAJ)
	
	# Open input file
	DefinedRegion=[]
	with open(f'{ResIDRegion}', 'r') as f:
		FileContent=f.readlines()
		
		for line in FileContent:
			if line.strip():
				DefinedRegion.append([int(line.split()[0]), int(line.split()[1]), line.split()[2]])


	# Calculate volume multiprocessing					
	pool = mp.Pool(mp.cpu_count())  
	pool.starmap(MeasuringVolume, [(u,ts,DefinedRegion,file_output) for ts in range(len(u.universe.trajectory))])
' > $WorkDir/RBS_volume_calculation/Volume_computations.py
	
# Create a .txt file that define interested region to compute volume
echo '91 91 A
129 134 A
141 141 A
149 149 A
151 152 A
179 182 A
186 186 A
189 191 A
218 218 A
220 226 A
' > $WorkDir/RBS_volume_calculation/InputFile.txt	
	
# Here: '91 91 A' indicates that all atoms of the amino acid id 91 belonging to chain A is part of the defined region	
# '129 134 A' indicates that all atoms of the amino acid id 129 - 134 belonging to chain A is also part of the defined region
	
# Measuring sampled RBS volume values over the course of simulations
python3 $WorkDir/RBS_volume_calculation/Volume_computations.py $WorkDir/ConcatenatedTraj/trajout.xtc $WorkDir/ConcatenatedTraj/frame0.pdb $WorkDir/RBS_volume_calculation/InputFile.txt $WorkDir/RBS_volume_calculation/VolumeStrctValues.txt

# Measuring RBS volume value of initial/starting structure
python3 $WorkDir/RBS_volume_calculation/Volume_computations.py $WorkDir/ConcatenatedTraj/frame0.pdb $WorkDir/ConcatenatedTraj/frame0.pdb $WorkDir/RBS_volume_calculation/InputFile.txt $WorkDir/RBS_volume_calculation/VolumeInitialStrctValue.txt
 

