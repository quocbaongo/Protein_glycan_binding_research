# Required input
WorkDir=/path/to/working/directory
ProteinHybridDir=/path/to/protein/hybrid/topology/stateA

TripeptideHybridDir=/path/to/protein/hybrid/topology/stateB

LIGStrct=/path/to/glycan/structure/in/pdb
LIGTop=/path/to/glycan/topology/in/itp
LIG=/glycan/name/in/glycan/topology/in/itp
GlycamFF=/path/to/glycam06h.itp
GMX=gmx_2021.2
PACKMOL=/path/to/packmol


		####################################################################################################################################



# Complex structure and topology preparation	
# Generate a complex of protein-receptor
mkdir $WorkDir/Complex_mut_hybrid_topol
cd $WorkDir/Complex_mut_hybrid_topol

# Generate complex of double topology protein and ligand in form of pdb file
grep -h ATOM $ProteinHybridDir/hybrid.pdb $LIGStrct >| $WorkDir/Complex_mut_hybrid_topol/Protein_LIG.pdb

# Renumber atom id of the pdb file of the complex structure
$GMX editconf -f $WorkDir/Complex_mut_hybrid_topol/Protein_LIG.pdb -o $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb

# Generate position restraint on the CA atom close to the center of geometry of the Protein_ligand complex
# Use python script to find the CA atom on the surface of protein that is close to the center of geometry of the Protein Ligand complex

echo "import MDAnalysis as mda
import numpy as np
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

if __name__ == '__main__':

	# Import PDB file
	pdb_file=sys.argv[1]

	# Load the PDB file using MDAnalysis
	u = mda.Universe(pdb_file)

	# Get protein and ligand selections 
	selection = 'protein or not protein'

	# Select the protein and ligand atoms
	Complex_center_geometry = u.select_atoms(selection).center_of_mass()
	
	# Protein CA atoms
	protein_ca_atoms = u.select_atoms('protein and name CA')
	
	# List of distance to center of geometry and
	# list of atom ID of CA atom of protein
	
	DistList=[]
	AtomIDList=[]

	
	for atom in protein_ca_atoms:
		DistList.append(MeasuringDistance(Complex_center_geometry,atom.position))
		AtomIDList.append(atom)
		
	
	ShortestDistIdx=DistList.index(min(DistList))
	print(AtomIDList[ShortestDistIdx].id)
" > $WorkDir/Complex_mut_hybrid_topol/Detect_CA_ID_COM.py

ProteinCAIDCOM=$(python3 $WorkDir/Complex_mut_hybrid_topol/Detect_CA_ID_COM.py $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb)


# Generate position restraint on the CA atom on the surface of protein 
#that is close to the center of geometry of the Protein Ligand complex

$GMX make_ndx -f $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb -o $WorkDir/Complex_mut_hybrid_topol/indexComplex.ndx <<EOF
a $ProteinCAIDCOM
q
EOF

echo "a_$ProteinCAIDCOM" | $GMX genrestr -f $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb -n $WorkDir/Complex_mut_hybrid_topol/indexComplex.ndx -o $WorkDir/Complex_mut_hybrid_topol/posre_COM_Protein.itp



echo "import MDAnalysis as mda
import numpy as np
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

if __name__ == '__main__':

	# Import PDB file
	pdb_file=sys.argv[1]

	# Load the PDB file using MDAnalysis
	u = mda.Universe(pdb_file)

	# Get protein and ligand selections 
	selection = 'protein and resid 1:320'

	# Select the protein and ligand atoms
	Complex_center_geometry = u.select_atoms(selection).center_of_mass()
	
	# Protein CA atoms
	protein_ca_atoms = u.select_atoms('protein and name CA')
	
	# List of distance to center of geometry and
	# list of atom ID of CA atom of protein
	
	DistList=[]
	AtomIDList=[]

	
	for atom in protein_ca_atoms:
		DistList.append(MeasuringDistance(Complex_center_geometry,atom.position))
		AtomIDList.append(atom)
		
	
	ShortestDistIdx=DistList.index(min(DistList))
	print(AtomIDList[ShortestDistIdx].id)
" > $WorkDir/Complex_mut_hybrid_topol/Detect_CA_ID_COM_2.py


ProteinCAIDCOM2=$(python3 $WorkDir/Complex_mut_hybrid_topol/Detect_CA_ID_COM_2.py $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb)


# Generate position restraint on the CA atom on the surface of protein 
#that is close to the center of geometry of the Protein Ligand complex

$GMX make_ndx -f $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb -n $WorkDir/Complex_mut_hybrid_topol/indexComplex.ndx -o $WorkDir/Complex_mut_hybrid_topol/indexComplex.ndx<<EOF
a $ProteinCAIDCOM2
"a_$ProteinCAIDCOM" | "a_$ProteinCAIDCOM2"
q
EOF

echo "a_$ProteinCAIDCOM2_a_$ProteinCAIDCOM" | $GMX genrestr -f $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb -n $WorkDir/Complex_mut_hybrid_topol/indexComplex.ndx -o $WorkDir/Complex_mut_hybrid_topol/posre_COM_Protein.itp -fc 1000 1000 1000



# Write to topology file
echo -e "\n#ifdef POSRES_TRANS\n#include \"$WorkDir/Complex_mut_hybrid_topol/posre_COM_Protein.itp\"\n#endif\n" >> $ProteinHybridDir/hybrid.itp




# Generate position restraint for ligand
$GMX make_ndx -f $LIGStrct -o $WorkDir/Complex_mut_hybrid_topol/index_LIG.ndx <<EOF
0 & ! a H*
q
EOF

echo 'System_&_!H*' | $GMX genrestr -f $LIGStrct -n $WorkDir/Complex_mut_hybrid_topol/index_LIG.ndx -o $WorkDir/Complex_mut_hybrid_topol/posre_LIG.itp -fc 1000 1000 1000

# Write position restraint to glycan topology file
echo -e "\n; Ligand position restraints\n#ifdef POSRES\n#include \"$WorkDir/Complex_mut_hybrid_topol/posre_LIG.itp\"\n#endif\n" >> $LIGTop


# Generate position restraint on the CA atom close to the center of geometry of the tripeptide
cd $TripeptideHybridDir

# Use python script to find the CA atom on the surface of protein that is close to the center of geometry of the Protein Ligand complex

echo "import MDAnalysis as mda
import numpy as np
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

if __name__ == '__main__':

	# Import PDB file
	pdb_file=sys.argv[1]

	# Load the PDB file using MDAnalysis
	u = mda.Universe(pdb_file)

	# Get protein and ligand selections 
	selection = 'protein or not protein'

	# Select the protein and ligand atoms
	Complex_center_geometry = u.select_atoms(selection).center_of_mass()
	
	# Protein CA atoms
	protein_ca_atoms = u.select_atoms('protein and name CA')
	
	# List of distance to center of geometry and
	# list of atom ID of CA atom of protein
	
	DistList=[]
	AtomIDList=[]

	
	for atom in protein_ca_atoms:
		DistList.append(MeasuringDistance(Complex_center_geometry,atom.position))
		AtomIDList.append(atom)
		
	
	ShortestDistIdx=DistList.index(min(DistList))
	print(AtomIDList[ShortestDistIdx].id)
" > $TripeptideHybridDir/Detect_CA_ID_COM.py

TripeptideCAIDAtom=$(python3 $TripeptideHybridDir/Detect_CA_ID_COM.py $TripeptideHybridDir/hybrid.pdb)
 
$GMX make_ndx -f $TripeptideHybridDir/hybrid.pdb -o $TripeptideHybridDir/indexTripeptide.ndx <<EOF
a $TripeptideCAIDAtom
q
EOF



echo "import MDAnalysis as mda
import numpy as np
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

if __name__ == '__main__':

	# Import PDB file
	pdb_file=sys.argv[1]

	# Load the PDB file using MDAnalysis
	u = mda.Universe(pdb_file)

	# Get protein and ligand selections 
	selection = 'protein and resid 1:320'

	# Select the protein and ligand atoms
	Complex_center_geometry = u.select_atoms(selection).center_of_mass()
	
	# Protein CA atoms
	protein_ca_atoms = u.select_atoms('protein and name CA')
	
	# List of distance to center of geometry and
	# list of atom ID of CA atom of protein
	
	DistList=[]
	AtomIDList=[]

	
	for atom in protein_ca_atoms:
		DistList.append(MeasuringDistance(Complex_center_geometry,atom.position))
		AtomIDList.append(atom)
		
	
	ShortestDistIdx=DistList.index(min(DistList))
	print(AtomIDList[ShortestDistIdx].id)
" > $TripeptideHybridDir/Detect_CA_ID_COM_2.py


TripeptideCAIDAtom2=$(python3 $TripeptideHybridDir/Detect_CA_ID_COM_2.py $TripeptideHybridDir/hybrid.pdb)


$GMX make_ndx -f $TripeptideHybridDir/hybrid.pdb -n $TripeptideHybridDir/indexTripeptide.ndx -o $TripeptideHybridDir/indexTripeptide.ndx <<EOF
a $ProteinCAIDCOM2
"a_$ProteinCAIDCOM" | "a_$ProteinCAIDCOM2"
q
EOF


echo "a_$ProteinCAIDCOM2_a_$ProteinCAIDCOM" | $GMX genrestr -f $TripeptideHybridDir/hybrid.pdb -n $TripeptideHybridDir/indexTripeptide.ndx -o $TripeptideHybridDir/posre_COM_tripeptide.itp


# Write position restraint to protein topology
echo -e "\n#ifdef POSRES_TRANS\n#include \"$TripeptideHybridDir/posre_COM_tripeptide.itp\"\n#endif\n" >> $TripeptideHybridDir/hybrid.itp



	####################################################################################################################################
# Preparing the system of protein-ligand complex and tripeptide

mkdir $WorkDir/2System1Box
mkdir $WorkDir/2System1Box/SystemPreparation
cd $WorkDir/2System1Box/SystemPreparation

# Prepare system of protein-ligand complex and tripeptide topology
cp $ProteinHybridDir/hybrid.top $WorkDir/2System1Box/SystemPreparation/system.top

cp $ProteinHybridDir/hybrid.itp $WorkDir/2System1Box/SystemPreparation
cp $ProteinHybridDir/hybrid_posre.itp $WorkDir/2System1Box/SystemPreparation

# Include the ligand topology to the system topology
cat $ProteinHybridDir/hybrid.top | sed "/forcefield\.itp\"/a\#include \"$GlycamFF\"" >| $WorkDir/2System1Box/SystemPreparation/system.top

sed -i "\|#include \"hybrid.itp\"|a\ \n#include \"$LIGTop\"\n" $WorkDir/2System1Box/SystemPreparation/system.top

echo "$LIG 1" >> $WorkDir/2System1Box/SystemPreparation/system.top


# Include the topology of the tripeptide to the system topology
sed -i "/#include \"amber99sb-star-ildn-mut.ff\/tip3p.itp\"/i #include \"$TripeptideHybridDir/hybrid.itp\"" $WorkDir/2System1Box/SystemPreparation/system.top

echo "Protein_chain_X 1" >> $WorkDir/2System1Box/SystemPreparation/system.top


# Prepare system of protein-ligand complex and tripeptide structure in form of pdb file
# Put protein-ligand complex in 1 simulation box with tripeptide using packmol
#Create complex_assembly.txt file (input for packmol)
echo "#minimum distance between protein-ligand complex and tripeptide
tolerance 40.0
#The file type of input and output files is PDB
filetype pdb
#The name of the output file
output $WorkDir/2System1Box/SystemPreparation/System.pdb
#distance from the edges of box
add_box_sides 1.0
structure $WorkDir/Complex_mut_hybrid_topol/Protein_LIG_renumbered.pdb 
 seed -1
 number 1
 inside box 0 0 0 120 120 120
end structure

structure $TripeptideHybridDir/hybrid.pdb
 seed -1
 number 1
 inside box 0 0 0 120 120 120
end structure
" > $WorkDir/2System1Box/SystemPreparation/complex.txt

$PACKMOL < $WorkDir/2System1Box/SystemPreparation/complex.txt




		####################################################################################################################################
# Start constructing simulation box
mkdir $WorkDir/2System1Box/SimulationBox
cd $WorkDir/2System1Box/SimulationBox

$GMX editconf -f $WorkDir/2System1Box/SystemPreparation/System.pdb -o $WorkDir/2System1Box/SimulationBox/System_box.pdb -bt dodecahedron -d 1.0 -c

$GMX solvate -cp $WorkDir/2System1Box/SimulationBox/System_box.pdb -cs spc216.gro -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/SimulationBox/System_water.pdb

$GMX grompp -f $WorkDir/MDPFile/genion.mdp -c $WorkDir/2System1Box/SimulationBox/System_water.pdb -r $WorkDir/2System1Box/SimulationBox/System_water.pdb -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/SimulationBox/tpr.tpr -maxwarn 1

echo "SOL" | $GMX genion -s $WorkDir/2System1Box/SimulationBox/tpr.tpr -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/SimulationBox/System_ions.pdb -conc 0.15 -pname NA -nname CL -neutral



# EM, Equilibration, Equilibrium production, Non-Equilibrium transition
#each process was repeated 3 times.
for i in 1 
do
	mkdir $WorkDir/2System1Box/Rep$i
done

# Energy minimization
# State A

for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/EM
	mkdir $WorkDir/2System1Box/Rep$i/EM/StateA
	cd $WorkDir/2System1Box/Rep$i/EM/StateA

	$GMX grompp -f $WorkDir/MDPFile/forward/equil_md/f_enmin.mdp -c $WorkDir/2System1Box/SimulationBox/System_ions.pdb -r $WorkDir/2System1Box/SimulationBox/System_ions.pdb -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/EM/StateA/tpr.tpr 
	
done

# State B
for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/EM/StateB
	cd $WorkDir/2System1Box/Rep$i/EM/StateB

	$GMX grompp -f $WorkDir/MDPFile/reverse/equil_md/r_enmin.mdp -c $WorkDir/2System1Box/SimulationBox/System_ions.pdb -r $WorkDir/2System1Box/SimulationBox/System_ions.pdb -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/EM/StateB/tpr.tpr 

done



for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i/EM/

	srun $GMX mdrun -s tpr.tpr -multidir State*
done



# NPT
# State A
for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/NPT
	mkdir $WorkDir/2System1Box/Rep$i/NPT/StateA
	cd $WorkDir/2System1Box/Rep$i/NPT/StateA
	

	$GMX grompp -f $WorkDir/MDPFile/forward/equil_md/f_npt.mdp -c $WorkDir/2System1Box/Rep$i/EM/StateA/confout.gro -r $WorkDir/2System1Box/Rep$i/EM/StateA/confout.gro -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/NPT/StateA/tpr.tpr -maxwarn 1

done


# State B
for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/NPT/StateB
	cd $WorkDir/2System1Box/Rep$i/NPT/StateB


	$GMX grompp -f $WorkDir/MDPFile/reverse/equil_md/r_npt.mdp -c $WorkDir/2System1Box/Rep$i/EM/StateB/confout.gro -r $WorkDir/2System1Box/Rep$i/EM/StateB/confout.gro -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/NPT/StateB/tpr.tpr -maxwarn 1

done 


cd $WorkDir/2System1Box

srun $GMX mdrun -s tpr.tpr -multidir Rep*/NPT/State*

# MD
# State A
for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/MD
	mkdir $WorkDir/2System1Box/Rep$i/MD/StateA
	cd $WorkDir/2System1Box/Rep$i/MD/StateA
	
	$GMX grompp -f $WorkDir/MDPFile/forward/equil_md/f_equil.mdp -c $WorkDir/2System1Box/Rep$i/NPT/StateA/confout.gro -r $WorkDir/2System1Box/Rep$i/NPT/StateA/confout.gro -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/MD/StateA/tpr.tpr -maxwarn 1

done

# State B
for i in 1
do 
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/MD/StateB
	cd $WorkDir/2System1Box/Rep$i/MD/StateB

	$GMX grompp -f $WorkDir/MDPFile/reverse/equil_md/r_equil.mdp -c $WorkDir/2System1Box/Rep$i/NPT/StateB/confout.gro -r $WorkDir/2System1Box/Rep$i/NPT/StateB/confout.gro -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/MD/StateB/tpr.tpr -maxwarn 1

done 



# MDRun
cd $WorkDir/2System1Box

srun $GMX mdrun -s tpr.tpr -multidir Rep*/MD/State*




		####################################################################################################################################	
	
# Non-Equilibrium transition simulation

# Extracting frames equidistantly along the equilibrium simulation trajectory 
# State A

for i in 1
do
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/FrameExtraction
	mkdir $WorkDir/2System1Box/Rep$i/FrameExtraction/StateA
	cd $WorkDir/2System1Box/Rep$i/FrameExtraction/StateA

	echo "System" | $GMX trjconv -f $WorkDir/2System1Box/Rep$i/MD/StateA/traj_comp.xtc -s $WorkDir/2System1Box/Rep$i/MD/StateA/tpr.tpr -o $WorkDir/2System1Box/Rep$i/FrameExtraction/StateA/frame_.pdb -ur compact -pbc mol -b 10000 -e 60000 -skip 5 -sep 

done


# State B

for i in 1
do
	cd $WorkDir/2System1Box/Rep$i
	mkdir $WorkDir/2System1Box/Rep$i/FrameExtraction/StateB
	cd $WorkDir/2System1Box/Rep$i/FrameExtraction/StateB

	echo "System" | $GMX trjconv -f $WorkDir/2System1Box/Rep$i/MD/StateB/traj_comp.xtc -s $WorkDir/2System1Box/Rep$i/MD/StateB/tpr.tpr -o $WorkDir/2System1Box/Rep$i/FrameExtraction/StateB/frame_.pdb -ur compact -pbc mol -b 10000 -e 60000 -skip 5 -sep 

done


# Prepare tpr files for the non-equilibrium transition btw WT->MT and MT->WT

for i in 1
do
	mkdir $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition
	mkdir $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateA
	mkdir $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateB	

done	

# State A
for i in 1
do
	for ((j=0;j<100;j++))
	do
		mkdir $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateA/frame$j
		cd $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateA/frame$j
		
		$GMX grompp -f $WorkDir/MDPFile/forward/nonequil_md/f_nonequil.mdp -c $WorkDir/2System1Box/Rep$i/FrameExtraction/StateA/frame_$j.pdb -r $WorkDir/2System1Box/Rep$i/FrameExtraction/StateA/frame_$j.pdb -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateA/frame$j/tpr.tpr -maxwarn 1

	done
done

# State B
for i in 1
do
	for ((j=0;j<100;j++))
	do
		mkdir $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateB/frame$j
		cd $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateB/frame$j

		$GMX grompp -f $WorkDir/MDPFile/reverse/nonequil_md/r_nonequil.mdp -c $WorkDir/2System1Box/Rep$i/FrameExtraction/StateB/frame_$j.pdb -r $WorkDir/2System1Box/Rep$i/FrameExtraction/StateB/frame_$j.pdb -p $WorkDir/2System1Box/SystemPreparation/system.top -o $WorkDir/2System1Box/Rep$i/NonEquilibriumTransition/StateB/frame$j/tpr.tpr -maxwarn 1
	
	done
done


# MD Run using parallelization scheme

