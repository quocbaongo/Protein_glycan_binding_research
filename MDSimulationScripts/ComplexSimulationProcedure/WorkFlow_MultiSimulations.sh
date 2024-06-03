# Required input
WorkDir=/path/to/working/directory
InputProteinStrct=/path/to/input/protein/structure/in/pdb
LIGStrct=/path/to/glycan/structure/in/pdb
LIGTop=/path/to/glycan/topology/in/itp
LIG=/glycan/name/in/glycan/topology/in/itp
GMX=gmx_2021.2


			############################ Initial MD simulation under high temperature ###################################



mkdir $WorkDir/HighEnergySampling
cd $WorkDir/HighEnergySampling

mkdir $WorkDir/HighEnergySampling/PDB2GMX
cd $WorkDir/HighEnergySampling/PDB2GMX

$GMX pdb2gmx -f $InputProteinStrct -o $WorkDir/HighEnergySampling/PDB2GMX/protein_pdb2gmx.pdb -p $WorkDir/HighEnergySampling/PDB2GMX/topol.top -ignh -merge all -ter <<EOF
8 
1
1
1
1
1
EOF

# The charmm36-jul2022 force field is selected
# Make N and C-terminal neutral backbone

# Create the structure and topology for protein-ligand complex
# Grep Protein LIG assembly
grep -h ATOM $WorkDir/HighEnergySampling/PDB2GMX/protein_pdb2gmx.pdb $LIGStrct >| $WorkDir/HighEnergySampling/PDB2GMX/Protein_Lig_complex.pdb

cp $WorkDir/HighEnergySampling/PDB2GMX/topol.top $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top

cat $WorkDir/HighEnergySampling/PDB2GMX/topol.top | sed "/forcefield\.itp\"/a\#include \"$LIGTop\"" >| $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top

echo "$LIG 1" >> $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top


# Generate Simulation box

mkdir $WorkDir/HighEnergySampling/Solvation
cd $WorkDir/HighEnergySampling/Solvation

$GMX editconf -f $WorkDir/HighEnergySampling/PDB2GMX/Protein_Lig_complex.pdb -o $WorkDir/HighEnergySampling/Solvation/ComplexBox.pdb -bt dodecahedron -d 1.0 -c

$GMX solvate -cp $WorkDir/HighEnergySampling/Solvation/ComplexBox.pdb -cs spc216.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -o $WorkDir/HighEnergySampling/Solvation/solv.gro

# Generate ions.mdp file
echo '; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = cutoff	; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
' > $WorkDir/HighEnergySampling/Solvation/ions.mdp

$GMX grompp -f $WorkDir/HighEnergySampling/Solvation/ions.mdp -c $WorkDir/HighEnergySampling/Solvation/solv.gro -r $WorkDir/HighEnergySampling/Solvation/solv.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -o $WorkDir/HighEnergySampling/Solvation/ions.tpr -maxwarn 1

echo "SOL" | $GMX genion -s $WorkDir/HighEnergySampling/Solvation/ions.tpr -o $WorkDir/HighEnergySampling/Solvation/solv_ions.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -pname NA -nname CL -neutral -conc 0.15

# Create LIG restraint
mkdir $WorkDir/HighEnergySampling/LIGRestraint
cd $WorkDir/HighEnergySampling/LIGRestraint

$GMX make_ndx -f $LIGStrct -o $WorkDir/HighEnergySampling/LIGRestraint/index_LIG.ndx <<EOF
0 & ! a H*
q
EOF

echo 'System_&_!H*' | $GMX genrestr -f $LIGStrct -n $WorkDir/HighEnergySampling/LIGRestraint/index_LIG.ndx -o $WorkDir/HighEnergySampling/LIGRestraint/posre_LIG.itp -fc 1000 1000 1000

# Edit $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top file
sed -i "\|#include \"$LIGTop\"|a\ \n; Ligand position restraints\n#ifdef POSRES\n#include \"$WorkDir/HighEnergySampling/LIGRestraint/posre_LIG.itp\"\n#endif" $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top


# EM
mkdir $WorkDir/HighEnergySampling/EM		
cd $WorkDir/HighEnergySampling/EM		
		
# Generate em.mdp file
echo "; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
" > $WorkDir/HighEnergySampling/EM/em.mdp

$GMX grompp -f $WorkDir/HighEnergySampling/EM/em.mdp -c $WorkDir/HighEnergySampling/Solvation/solv_ions.gro -r $WorkDir/HighEnergySampling/Solvation/solv_ions.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -o $WorkDir/HighEnergySampling/EM/em.tpr
		
srun $GMX mdrun -s $WorkDir/HighEnergySampling/EM/em.tpr 		


# NVT
mkdir $WorkDir/HighEnergySampling/NVT
cd $WorkDir/HighEnergySampling/NVT		

# Make index file for tc-grps
# NVT at 398K for the system to overcome energy barrier

$GMX make_ndx -f $WorkDir/HighEnergySampling/EM/em.tpr -o $WorkDir/HighEnergySampling/NVT/index.ndx <<EOF
q
EOF


echo "import MDAnalysis as mda
import sys

if __name__ == '__main__':

	# Input file	
	PDBFile=sys.argv[1]

	# Load the PDB file using MDAnalysis
	u = mda.Universe(PDBFile)

	# Access the AtomGroup representing all atoms in the PDB file
	atoms = u.atoms	
	
	# Extract the atom id
	AtomIDList=[]
	for atom in atoms:
		AtomIDList.append(atom.id)

	# Write atom id to index file
	out_file = open('index.ndx', 'a')
	out_file.write('[ Protein_LIG ]\n')
					
	for i in range(0, len(AtomIDList), 15):
		for idex in AtomIDList[i:i+15]:
			if idex < 10:
				out_file.write('   '+str(idex) + ' ')
			
			elif (idex > 9) and (idex < 100):
				out_file.write('  '+str(idex) + ' ')
				
			elif (idex > 99) and (idex < 1000):
				out_file.write(' '+str(idex) + ' ')
			else:
				out_file.write(''+str(idex) + ' ')										
				
		out_file.write('\n')	
" > $WorkDir/HighEnergySampling/NVT/IndexGenerator.py

python3 $WorkDir/HighEnergySampling/NVT/IndexGenerator.py $WorkDir/HighEnergySampling/Solvation/ComplexBox.pdb

# Generate nvt.mdp file
echo "title                   = NVT equilibration 
define                  = -DPOSRES  ; position restrain 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 250000     ; 2 * 250000 = 500 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500   ; save energies every 1.0 ps
nstlog                  = 500   ; update log file every 1.0 ps
nstxout-compressed      = 500   ; save coordinates every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_LIG Water_and_ions  
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 398   398                            ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 398       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
" > $WorkDir/HighEnergySampling/NVT/nvt.mdp


$GMX grompp -f $WorkDir/HighEnergySampling/NVT/nvt.mdp -c $WorkDir/HighEnergySampling/EM/confout.gro -r $WorkDir/HighEnergySampling/EM/confout.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -n $WorkDir/HighEnergySampling/NVT/index.ndx -o $WorkDir/HighEnergySampling/NVT/nvt.tpr

srun $GMX mdrun -s $WorkDir/HighEnergySampling/NVT/nvt.tpr 


# NPT
mkdir $WorkDir/HighEnergySampling/NPT
cd $WorkDir/HighEnergySampling/NPT

cp $WorkDir/HighEnergySampling/NVT/index.ndx $WorkDir/HighEnergySampling/NPT

# Generate npt.mdp file
echo "title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 500000    ; 2 * 500000 = 1000 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_LIG Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 398   398                            ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = Berendsen                     ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT
" > $WorkDir/HighEnergySampling/NPT/npt.mdp


$GMX grompp -f $WorkDir/HighEnergySampling/NPT/npt.mdp -c $WorkDir/HighEnergySampling/NVT/confout.gro -r $WorkDir/HighEnergySampling/NVT/confout.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -n $WorkDir/HighEnergySampling/NPT/index.ndx -o $WorkDir/HighEnergySampling/NPT/npt.tpr

srun $GMX mdrun -s $WorkDir/HighEnergySampling/NPT/npt.tpr 



# Simulate system at 398K for 1ns
# MD
mkdir $WorkDir/HighEnergySampling/MD
cd $WorkDir/HighEnergySampling/MD

cp $WorkDir/HighEnergySampling/NPT/index.ndx $WorkDir/HighEnergySampling/MD

# Generate md.mdp file
echo "title                   = Protein-ligand complex MD simulation 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 500000    ; 2 * 500000 = 1000 ps (1 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save coordinates every 10.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_LIG Water_and_ions    ; two coupling groups - more accurate

tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 398   398                            ; reference temperature, one for each group, in K

; Pressure coupling 
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration 

" > $WorkDir/HighEnergySampling/MD/md.mdp


$GMX grompp -f $WorkDir/HighEnergySampling/MD/md.mdp -c $WorkDir/HighEnergySampling/NPT/confout.gro -r $WorkDir/HighEnergySampling/NPT/confout.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -n $WorkDir/HighEnergySampling/MD/index.ndx -o $WorkDir/HighEnergySampling/MD/md.tpr

srun $GMX mdrun -s $WorkDir/HighEnergySampling/MD/md.tpr 



# Fix the PBC
mkdir $WorkDir/HighEnergySampling/MD/FixPBC
cd $WorkDir/HighEnergySampling/MD/FixPBC

echo "Protein_LIG" "Protein_LIG" | $GMX trjconv -f $WorkDir/HighEnergySampling/MD/traj_comp.xtc -s $WorkDir/HighEnergySampling/MD/md.tpr -o $WorkDir/HighEnergySampling/MD/FixPBC/md_center.xtc -center -pbc mol -ur compact -n $WorkDir/HighEnergySampling/MD/index.ndx

echo "Protein_LIG" "Protein_LIG" | $GMX trjconv -f $WorkDir/HighEnergySampling/MD/FixPBC/md_center.xtc -s $WorkDir/HighEnergySampling/MD/md.tpr -o $WorkDir/HighEnergySampling/MD/FixPBC/md_center_fit.xtc -fit rot+trans -n $WorkDir/HighEnergySampling/MD/index.ndx

echo "Protein_LIG" | $GMX trjconv -f $WorkDir/HighEnergySampling/MD/FixPBC/md_center_fit.xtc -s $WorkDir/HighEnergySampling/MD/md.tpr -dump 0 -o $WorkDir/HighEnergySampling/MD/FixPBC/frame0.pdb -n $WorkDir/HighEnergySampling/MD/index.ndx

rm $WorkDir/HighEnergySampling/MD/FixPBC/md_center.xtc


				#########################################################################################################

# From the the complex simulation at high temperature
# 10 conformations were extracted equidistantly
# which subsequently would serve as the starting conformations
# for the multiple independent simulations

mkdir $WorkDir/StartingConformations
cd $WorkDir/StartingConformations

echo "System" | $GMX trjconv -f $WorkDir/HighEnergySampling/MD/traj_comp.xtc -s $WorkDir/HighEnergySampling/MD/md.tpr -o $WorkDir/StartingConformations/frame_.pdb -sep -skip 11 -ur compact -pbc mol



# Sampling simulation
mkdir $WorkDir/MultiSimulations
cd $WorkDir/MultiSimulations

#Create directory for each initial conformation
for((i=0;i<10;i++))
do
	mkdir $WorkDir/MultiSimulations/frame_$i
done


# Generate .mdp file for multiple independent simulations
mkdir $WorkDir/MultiSimulations/MDPFile
cd $WorkDir/MultiSimulations/MDPFile

# Generate nvt.mdp file
echo "title                   = NVT equilibration 
define                  = -DPOSRES  ; position restrain 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 250000     ; 2 * 250000 = 500 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500   ; save energies every 1.0 ps
nstlog                  = 500   ; update log file every 1.0 ps
nstxout-compressed      = 500   ; save coordinates every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_LIG Water_and_ions    ; two coupling groups - more accurate

tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                            ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
" > $WorkDir/MultiSimulations/MDPFile/nvt.mdp

# Generate npt.mdp file
echo "title                   = NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 500000    ; 2 * 500000 = 1000 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_LIG Water_and_ions    ; two coupling groups - more accurate

tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                            ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = Berendsen                     ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT
" > $WorkDir/MultiSimulations/MDPFile/npt.mdp



# Generate md.mdp file
echo "title                   = Protein-ligand complex MD simulation 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000000  ; 2 * 50000000 = 100000 ps (100 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save coordinates every 10.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_LIG Water_and_ions    ; two coupling groups - more accurate

tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                            ; reference temperature, one for each group, in K

; Pressure coupling 
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration 
" > $WorkDir/MultiSimulations/MDPFile/md.mdp


# NVT and assign initial velocities
for((i=0;i<10;i++))
do
	mkdir $WorkDir/MultiSimulations/frame_$i/NVT
	cd $WorkDir/MultiSimulations/frame_$i/NVT

	cp $WorkDir/HighEnergySampling/NVT/index.ndx $WorkDir/MultiSimulations/frame_$i/NVT

	$GMX grompp -f $WorkDir/MultiSimulations/MDPFile/nvt.mdp -c $WorkDir/StartingConformations/frame_$i.pdb -r $WorkDir/StartingConformations/frame_$i.pdb -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -n $WorkDir/MultiSimulations/frame_$i/NVT/index.ndx -o $WorkDir/MultiSimulations/frame_$i/NVT/nvt.tpr

done


# Srun
cd $WorkDir/MultiSimulations

srun $GMX mdrun -s nvt.tpr -multidir $WorkDir/MultiSimulations/frame_*/NVT



# NPT

for((i=0;i<10;i++))
do
	mkdir $WorkDir/MultiSimulations/frame_$i/NPT
	cd $WorkDir/MultiSimulations/frame_$i/NPT

	cp $WorkDir/HighEnergySampling/NVT/index.ndx $WorkDir/MultiSimulations/frame_$i/NPT

	$GMX grompp -f $WorkDir/MultiSimulations/MDPFile/npt.mdp -c $WorkDir/MultiSimulations/frame_$i/NVT/confout.gro -r $WorkDir/MultiSimulations/frame_$i/NVT/confout.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -n $WorkDir/MultiSimulations/frame_$i/NPT/index.ndx -o $WorkDir/MultiSimulations/frame_$i/NPT/npt.tpr
	
done

# Srun
srun $GMX mdrun -s npt.tpr -multidir $WorkDir/MultiSimulations/frame_*/NPT


# MD run

for((i=0;i<10;i++))
do
	mkdir $WorkDir/MultiSimulations/frame_$i/MD
	cd $WorkDir/MultiSimulations/frame_$i/MD
	
	cp $WorkDir/HighEnergySampling/NVT/index.ndx $WorkDir/MultiSimulations/frame_$i/MD

	$GMX grompp -f $WorkDir/MultiSimulations/MDPFile/md.mdp -c $WorkDir/MultiSimulations/frame_$i/NPT/confout.gro -r $WorkDir/MultiSimulations/frame_$i/NPT/confout.gro -p $WorkDir/HighEnergySampling/PDB2GMX/Protein_LIG.top -n $WorkDir/MultiSimulations/frame_$i/MD/index.ndx -o $WorkDir/MultiSimulations/frame_$i/MD/md.tpr

done

# Production run
srun $GMX mdrun -s md.tpr -multidir $WorkDir/MultiSimulations/frame_*/MD


