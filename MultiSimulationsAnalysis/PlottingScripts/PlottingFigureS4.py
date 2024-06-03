import os
import json
import sys
import argparse
import MDAnalysis
import matplotlib.pyplot as plt

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

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

if __name__ == "__main__":

    # flag list
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description = "Amino acids that are actively involved in\
                                                H-bond interactions with a sugar of glycan in the\
                                                researched complex")
                                                        
    parser.add_argument("--Input", type=validate_file, help="Input json file for number of hbond\
                                                        contacts with sugar for each amino acid\
                                                        (HydrogenBondCountDict.json)")
    parser.add_argument("--PDB", type=validate_file, help="Complex PDB file")
    parser.add_argument("--Output", action="store", help="Name of output png file (Optional)")
    parser.add_argument("--y_lower_limit", action="store", help="Lower limit of y axis (Optional)")
    parser.add_argument("--y_upper_limit", action="store", help="Upper limit of y axis (Optional)")    

    args = parser.parse_args()
    # Input file
    Input=args.Input   
    PDB=args.PDB

    # Load topology and trajectory to MDAnalysis
    uPDB = MDAnalysis.Universe(f"{PDB}")
    
    Database={}
    
    for resid,resname in zip(uPDB.select_atoms("protein").resids,uPDB.select_atoms("protein").resnames):
    	Database.update({resid:d[resname]})
    	

    # Open json file
    FileInput = open(Input)
    FileInputContent = json.load(FileInput)
   
    # Normalized dict
    TotalHbonds=0
    for key, value in FileInputContent.items():
        TotalHbonds+=FileInputContent[key]
        
    for key, value in FileInputContent.items():
        FileInputContent[key] /= TotalHbonds
        
    ########################################## Start Plotting ######################################################
    # Plotting
    # Define size of figure
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)

    # x axis ticks
    xticks=[]
    
    for key in FileInputContent.keys():
    	if abs(FileInputContent[key]) > 0.01:
    		xticks.append(str(key))
    		
    xticksLabel=[mem for mem in xticks]
    for i in range(len(xticksLabel)):
        xticksLabel[i] = f"{Database[int(xticksLabel[i])]}{xticksLabel[i]}"
        
    ax.bar(xticksLabel,[float(FileInputContent[key]) for key in xticks])
    
    # Value above each bar
    for p in ax.patches:
    	b = p.get_bbox()
    	if b.y1 < 0:
    		BarVal="{:+.3f}".format(round(b.y1,4))
    		ax.annotate(BarVal,((b.x0 + b.x1)/2,0),rotation=60,size=6,weight='bold')
    	else:
    		BarVal="{:+.3f}".format(round(b.y1,4))
    		ax.annotate(BarVal,((b.x0 + b.x1)/2,b.y1),rotation=60,size=6,weight='bold')    
    
    
    ax.set_xticks(xticksLabel, rotation=30)
    
    # x,yticks
    if args.y_lower_limit and args.y_upper_limit:
    	ax.set_ylim((float(args.y_lower_limit), float(args.y_upper_limit)))
    else:
    	pass
    	
    # Set label
    ax.set_xlabel("Residue ID",labelpad=10)
    ax.set_ylabel("H-bond Formation \nprobability",labelpad=10)

    plt.tight_layout()
    
    if args.Output:
    	plt.savefig(f"{args.Output}.png",dpi=600,bbox_inches="tight")
    else:
    	plt.savefig("Out.png",dpi=600,bbox_inches="tight")
    
