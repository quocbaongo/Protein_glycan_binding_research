import matplotlib.pyplot as plt
import argparse
import json
import os

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
    parser = argparse.ArgumentParser(description = "Detecting distance between center of geometry of \
                                                a glycan and CA atom of a target protein residue \
                                                in each time frame of the simulation trajectory")
                                                        
    parser.add_argument("--Input", type=validate_file, help="Input json file that contains the distance between center of geometry of \
    								a glycan and CA atom of a target protein residue \
    								in each time frame of the simulation trajectory")
    parser.add_argument("--y_lower_limit", action="store", help="Lower limit of y axis (Optional)")
    parser.add_argument("--y_upper_limit", action="store", help="Upper limit of y axis (Optional)")
    parser.add_argument("--Output", action="store", help="Name of output png file (Optional)")
    
    args = parser.parse_args()
    
    # Input file
    InFile=args.Input
        								
    # Open json file
    FileData=open(InFile)
    FileDataContent=json.load(FileData)
    FileDataContentSorted=sorted(FileDataContent.items(), key=lambda t: get_key(t[0]))
    FileDataContentSorted=[(str(float(element[0])/100),element[1]/10.0) for element in FileDataContentSorted]
    
    
    # Prepare xticks
    xticks=[]
    for element in [i[0] for i in FileDataContentSorted]:
    	if float(element) % 10.0 == 0.0:
    		xticks.append(element)
    	else:
    		xticks.append('')    
 
 
    # Plotting Evolution
    plt.figure(figsize=(8,4))
    plt.plot([i[0] for i in FileDataContentSorted],[i[1] for i in FileDataContentSorted])
    plt.xticks(xticks)
    
    if args.y_lower_limit and args.y_upper_limit:
    	plt.ylim((float(args.y_lower_limit), float(args.y_upper_limit)))
    else:
    	pass
    	
    plt.xlabel('Time (ns)',size=10)
    plt.ylabel('Distance (nm)',size=10)
    
    if args.Output:
    	plt.savefig(f"{args.Output}.png",dpi=600,bbox_inches="tight")
    else:
    	plt.savefig("Out.png",dpi=600,bbox_inches="tight")    
    	
    								

