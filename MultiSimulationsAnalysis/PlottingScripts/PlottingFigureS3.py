import os
import json
import sys
import argparse
import numpy as np
import MDAnalysis
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
    return arg.split(",")

if __name__ == "__main__":

    # flag list
    parser = argparse.ArgumentParser()
    
    # Add an argument for the list of strings
    parser.add_argument("--InputEnergyFiles", metavar="N", nargs="+", help="Input Energy at each time frame files")
    parser.add_argument("--InputDistanceFiles", metavar="N", nargs="+", help="Input center of mass distance at each time frame files")
    parser.add_argument("--Output", action="store", help="Name of output png file (Optional)")

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Extracting distances
    DistData=[]
    for file in args.InputDistanceFiles:
        FileData=open(f"{file}")
        FileDataContent=json.load(FileData)
        
        for key,value in FileDataContent.items():
            DistData.append(float(value))
    
    DistData=[element/10.0 for element in DistData]
    
    # Extracting Coulombic and LJ energies
    LJData=[]
    ColData=[]
    
    for file in args.InputEnergyFiles:
        FileData=open(f"{file}")
        
        Temp=[]
        for line in FileData:
            if not (line.startswith("#") or line.startswith("@")):
                ColData.append(float(line.split()[1]))
                LJData.append(float(line.split()[2]))
    
    ########################################## Start Plotting ######################################################
    # Plotting
    # Define size of figure
    fig = plt.figure(figsize=(14, 7))
    gs = gridspec.GridSpec(20, 46)
    # Define the positions of the subplots.
    ax0=fig.add_subplot(gs[:15, :20])
    ax1=fig.add_subplot(gs[:15:, 26:])
        
    # Plotting
    # LJ vs distance plotting
    ax0.scatter(DistData, LJData)    
    
    # Coulomb vs distance plotting    
    ax1.scatter(DistData, ColData)
    
    # Set x axis ticks range
    ax0.set_xticks(np.arange(0.0, max(DistData)+1, 2.5))
    ax1.set_xticks(np.arange(0.0, max(DistData)+1, 2.5))
    ax0.set_xticklabels(ax0.get_xticks(),rotation=30)
    ax1.set_xticklabels(ax1.get_xticks(),rotation=30)
    
    # Set label
    ax0.set_xlabel("Distance (nm)",labelpad=10)
    ax1.set_xlabel("Distance (nm)",labelpad=10)
    ax0.set_ylabel("LJ Energy (kJ/mol)",labelpad=10)
    ax1.set_ylabel("Coulomb Energy (kJ/mol)",labelpad=10)
    
    # Set fig label
    ax0.set_title("A", loc="left", fontdict={"fontsize": 20})
    ax1.set_title("B", loc="left", fontdict={"fontsize": 20})
    
    plt.tight_layout()
    if args.Output:
        plt.savefig(f"{args.Output}.png",dpi=600,bbox_inches="tight")
    else:
        plt.savefig("Out.png",dpi=600,bbox_inches="tight")







