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

# Define a custom argument type for a list of strings
def list_of_strings(arg):
    return arg.split(',')

if __name__ == "__main__":

    # flag list
    parser = argparse.ArgumentParser()
    
    # Add an argument for the list of strings
    parser.add_argument("--InputRMSD", type=validate_file, help="Input RMSDs")
    parser.add_argument("--InitialRMSD", action="store", metavar="N", nargs="+", help="List of RMSD files for starting structures (Optional)")
    parser.add_argument("--InputAngle", type=validate_file, help="Input Angles")
    parser.add_argument("--InitialAngle", action="store", metavar="N", nargs="+", help="List of Angle files for starting structures (Optional)")
    parser.add_argument("--rmsd_lower_limit", action="store", help="Lower limit of RMSD (Optional)")
    parser.add_argument("--rmsd_upper_limit", action="store", help="Upper limit of RMSD (Optional)")
    parser.add_argument("--angle_lower_limit", action="store", help="Lower limit of angle (Optional)")
    parser.add_argument("--angle_upper_limit", action="store", help="Upper limit of angle (Optional)")
    parser.add_argument("--Output", action="store", help="Name of output png file (Optional)")


    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Open Input RMSD
    with open(f"{args.InputRMSD}", "r") as f:
        InputRMSDContent = f.readlines()

    RMSDList=[]
    for line in InputRMSDContent:
        if not (line.startswith('#') or line.startswith('@')):
            temp=[float(line.split()[0]), float(line.split()[1])]
            RMSDList.append(temp)
        
    SortedRMSDList=sorted(RMSDList, key=lambda x: x[0])

    # Finding RMSD original value
    InitialRMSD={}
    for file in args.InitialRMSD:
        if file == args.InputRMSD:
            continue
        else:
            # Extract name for key
            RMSDKey=""
            for i in os.path.basename(file).split():
                for letter in i:   
                    if letter.isdigit():
                        RMSDKey+=letter

            with open(f"{file}", 'r') as f:
                FileContent = f.readlines()
            for line in FileContent:
                if not (line.startswith('#') or line.startswith('@')):
                    InitialRMSDValue=float(line.split()[1])
            
            
            InitialRMSD.update({int(RMSDKey):InitialRMSDValue})
    

    # Open Input Angle
    with open(f"{args.InputAngle}", "r") as f:
        InputAngleContent = f.readlines()

    AngleList=[]
    for line in InputAngleContent:
        if not (line.startswith('#') or line.startswith('@')):
            temp=[float(line.split()[0]), float(line.split()[1])]
            AngleList.append(temp)
        
    SortedAngleList=sorted(AngleList, key=lambda x: x[0])
    
    # Finding angle original value
    InitialAngle={}
    for file in args.InitialAngle:
        if file == args.InputAngle:
            continue
        else:
            # Extract name for key
            AngleKey=""
            for i in os.path.basename(file).split():
                for letter in i:   
                    if letter.isdigit():
                        AngleKey+=letter
                    
            #print(AngleKey)
        
            with open(f"{file}", 'r') as f:
                FileContent = f.readlines()
            for line in FileContent:
                if not (line.startswith('#') or line.startswith('@')):
                    InitialAngleValue=float(line.split()[1])
            
            
            InitialAngle.update({int(AngleKey):InitialAngleValue})
    
    ########################################## Start Plotting ######################################################
    plt.figure(figsize=(8,4))
    
    # RMSD limit
    if args.rmsd_lower_limit and args.rmsd_upper_limit:
        RMSDLowLimit=float(args.rmsd_lower_limit)
        RMSDUpLimit=float(args.rmsd_upper_limit)
    else:
    	pass
    	
    # Angle limit
    if args.angle_lower_limit and args.angle_upper_limit:
        AngleLowLimit=float(args.angle_lower_limit)
        AngleUpLimit=float(args.angle_upper_limit)
    else:
    	pass    

    # Plotting 
    if ("RMSDLowLimit" and "RMSDUpLimit" and "AngleLowLimit" and "AngleUpLimit") in globals():
        plt.hist2d([i[1] for i in SortedRMSDList], [i[1] for i in SortedAngleList], bins=100, cmap='hot',range=[[RMSDLowLimit,RMSDUpLimit],[AngleLowLimit,AngleUpLimit]])
    else:
        plt.hist2d([i[1] for i in SortedRMSDList], [i[1] for i in SortedAngleList], bins=100, cmap='hot')
            
    # Initial positions:
    for key in InitialRMSD.keys():
        if key == min(list(InitialRMSD.keys())):
            plt.plot(InitialRMSD[key],InitialAngle[key],marker="o",color="white")
        else:
            plt.plot(InitialRMSD[key],InitialAngle[key],marker="o",color="blue")
	        
	# Set label
    plt.xlabel('RMSD (nm)',labelpad=10)
    plt.ylabel('Angle θ (degree)',labelpad=10)
    
    plt.tight_layout()
    if args.Output:
        plt.savefig(f"{args.Output}.png",dpi=600,bbox_inches="tight")
    else:
        plt.savefig("Out.png",dpi=600,bbox_inches="tight")       
    
