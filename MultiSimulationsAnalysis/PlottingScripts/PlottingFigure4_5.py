import os
import json
import sys
import argparse
import MDAnalysis
import matplotlib.pyplot as plt

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
    parser = argparse.ArgumentParser(description = "Comparison of hydrogen bond networks involving the target amino acid \
    							between the wild-type and mutant complexes")
                                                        
    parser.add_argument("--WTInput", type=validate_file, help="WT Input json file (HydrogenBondCountDict.json)")
    parser.add_argument("--WTTotalFrame", help="Total number of frame in WT trajectory, normalization purpose")
    parser.add_argument("--MTInput", type=validate_file, help="MT Input json file (HydrogenBondCountDict.json)")
    parser.add_argument("--MTTotalFrame", help="Total number of frame in MT trajectory, normalization purpose")
    parser.add_argument("--WTPDB", type=validate_file, help="WT PDB file")
    parser.add_argument("--Output", action="store", help="Name of output png file (Optional)")
    parser.add_argument("--y_lower_limit", action="store", help="Lower limit of y axis (Optional)")
    parser.add_argument("--y_upper_limit", action="store", help="Upper limit of y axis (Optional)")    

    args = parser.parse_args()
    # Input file
    WTInput=args.WTInput
    WTTotalFrame=args.WTTotalFrame    
    MTInput=args.MTInput
    MTTotalFrame=args.MTTotalFrame
    WTPDB=args.WTPDB


    # Load topology and trajectory to MDAnalysis
    uPDB_WT = MDAnalysis.Universe(f"{WTPDB}")
    
    WTDatabase={}
    
    
    for resid,resname,chainID in zip(uPDB_WT.select_atoms("all").resids,uPDB_WT.select_atoms("all").resnames,uPDB_WT.select_atoms("all").segids):
    	WTDatabase.update({f"{resid}{chainID}":resname})
    	
    # Open json file
    FileWTInput = open(WTInput)
    FileWTInputContent = json.load(FileWTInput)

    FileMTInput = open(MTInput)
    FileMTInputContent = json.load(FileMTInput)    
    # Normalized dict
    # WT
    for key, value in FileWTInputContent.items():
    	FileWTInputContent[key] /= int(WTTotalFrame)
    	
    # MT
    for key, value in FileMTInputContent.items():
    	FileMTInputContent[key] /= int(MTTotalFrame)
    	
    # Convert list to dict
    DiffDict={}

    # WT/MT HA23 Res94
    for key, value in FileWTInputContent.items():
    	Diff=FileWTInputContent[key]-FileMTInputContent[key]
    	DiffDict.update({key:Diff})
    

    ########################################## Start Plotting ######################################################
    # Plotting
    # Define size of figure
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    
    # x axis ticks
    xticks=[]
    
    for key in DiffDict.keys():
    	if abs(DiffDict[key]) >= 0.001:
    		xticks.append(str(key))
    		
    xticksLabel=[mem for mem in xticks]
    
    for i in range(len(xticksLabel)):	
    	xticksLabel[i] = f"{WTDatabase[xticksLabel[i]]}{xticksLabel[i]}"
    	

    ax.bar(xticksLabel,[float(DiffDict[key]) for key in xticks])
    
    # Value above each bar
    for p in ax.patches:
    	b = p.get_bbox()
    	if b.y1 < 0:
    		BarVal="{:+.3f}".format(round(b.y1,4))
    		ax.annotate(BarVal,((b.x0 + b.x1)/2,0),rotation=60,size=6,weight='bold')
    	else:
    		BarVal="{:+.3f}".format(round(b.y1,4))
    		ax.annotate(BarVal,((b.x0 + b.x1)/2,b.y1),rotation=60,size=6,weight='bold')    
    
    ax.set_xticklabels(xticksLabel, rotation=30)
    
    # x,yticks
    if args.y_lower_limit and args.y_upper_limit:
    	ax.set_ylim((float(args.y_lower_limit), float(args.y_upper_limit)))
    else:
    	pass
    	
    # Set label
    ax.set_xlabel("Residue ID",labelpad=10)
    ax.set_ylabel("H-bond Formation Propensity \n difference",labelpad=10)

    plt.tight_layout()
    
    if args.Output:
    	plt.savefig(f"{args.Output}.png",dpi=600,bbox_inches="tight")
    else:
    	plt.savefig("Out.png",dpi=600,bbox_inches="tight")
    	














    
    
    
    
    
    
