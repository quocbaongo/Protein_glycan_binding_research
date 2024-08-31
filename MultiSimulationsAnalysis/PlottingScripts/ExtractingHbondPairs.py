import os
import json
import argparse
import MDAnalysis
import numpy as np
from collections import Counter

def validate_file(f):

 	if not os.path.exists(f):
 		# Argparse uses the ArgumentTypeError to give a rejection message like:
 		# error: argument input: x does not exist
 		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
 
 	return f
    
def get_key(key):
 	try:
 		return int(key)
 	except ValueError:
 		return key
		
AminoAcid=["VAL","ILE","LEU","GLU","GLN","ASP",
	"ASN","HIS","TRP","PHE","TYR","ARG","LYS",
	"SER","THR","MET","ALA","GLY","PRO","CYS"]


if __name__ == "__main__":

	# flag list
	parser = argparse.ArgumentParser(description = "Extracting hydrogen bond pairs and probability of occuring from file containing detected hydrogen bonds per time step")
    
	# Add an argument for the list of strings
	parser.add_argument("--InputFile", type=validate_file, help="File containing detected hydrogen bonds per time step (often named HydrogenBondPerTimeStep.txt)", required=True)
	parser.add_argument("--PDBFile", type=validate_file, help="Structure file in .pdb format", required=True)
	parser.add_argument("--resid1", help="First input residue ID (residue id based on numbering of input structure defined in flag --PDBFile)", required=True)
	parser.add_argument("--chain1", help="Chain ID of first input residue ID defined in flag --resid1 (chain id based on chain numbering of input structure defined in flag --PDBFile)", required=True)
	parser.add_argument("--resid2", help="Second input residue ID (residue id based on numbering of input structure defined in flag --PDBFile)", required=True)
	parser.add_argument("--chain2", help="Chain ID of second input residue ID defined in flag --resid1 (chain id based on chain numbering of input structure defined in flag --PDBFile)", required=True)
	parser.add_argument("--Output", action="store", help="Name of output file. The output contains pair", required=True)
	
	# Parse the command-line arguments
	args = parser.parse_args()
 	
	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	uPDB = MDAnalysis.Universe(args.PDBFile)

	# Create dictionary for counting hydrogen bond
	HydrogenBondPerProteinRes={}
	for ProteinRes in np.unique([atom.resid for atom in uPDB.select_atoms("protein")]):
		HydrogenBondPerProteinRes.update({str(ProteinRes):0})
 		
	# Read HydrogenFile
	with open(args.InputFile) as f:
		HydrogenFileContent = f.read().splitlines()
 		
	Interaction=[]	
	# Read each hydrogen bond
	for line in HydrogenFileContent:
		HydrogenBondConnection=[]
		for atom in line.split():
			HydrogenBondConnection.append(atom.split("_")[0])
				
		HydrogenBondConnection=np.unique(HydrogenBondConnection)
		
		if (f"{args.resid1}{args.chain1}" in HydrogenBondConnection) and (f"{args.resid2}{args.chain2}" in HydrogenBondConnection):
			Interaction.append(line)


	CountDict=dict(Counter(Interaction))
	
	# Hbond probability
	TotalHbond=0
	
	for key,value in CountDict.items():
		TotalHbond+=value

	# Convert to percentage
	for key,value in CountDict.items():
		CountDict[key] /= TotalHbond

	# Write to file
	with open(f"{args.Output}.txt", "a") as f:
		for key,value in CountDict.items():
			f.write(f"{key} {value}\n")



