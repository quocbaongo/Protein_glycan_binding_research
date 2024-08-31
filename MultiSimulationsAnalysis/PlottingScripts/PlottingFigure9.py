(insilico_methodologies_2022) [ngoquoc1@mahti-login14 RBS_volume_plotting]$ cat VolumeValueDistribution.py 
import numpy as np
import matplotlib.pyplot as plt 	
import sys
import argparse

if __name__ == '__main__':

	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Plotting sampled volume values distribution')
	parser.add_argument('--WTFile', help='Sampled volume of WT structure')
	parser.add_argument('--WTOriFile', help='Original volume of WT structure')	
	parser.add_argument('--MTFile', help='Sampled volume of MT structure')
	parser.add_argument('--MTOriFile', help='Original volume of MT structure')
	parser.add_argument('--file_out', help='Name of output image in .png and .eps format')

	args = parser.parse_args()
	fileWT=args.WTFile
	fileMT=args.MTFile
	FileOriginalWT=args.WTOriFile	
	FileOriginalMT=args.MTOriFile
	FileOut=args.file_out

	# Open JsonFile for initial volume value
	WTInitialValue=0.0
	MTInitialValue=0.0	
	
	# WT 
	with open(FileOriginalWT) as text_file:
		for line in text_file:
			WTInitialValue=float(line.split()[1])
			
	# MT 
	with open(FileOriginalMT) as text_file:
		for line in text_file:
			MTInitialValue=float(line.split()[1])
	
	# Handling volume value in structural ensemble
	#WTHandling
	fileWT = open(fileWT, "r")
	WTTrajListConcat = fileWT.readlines()
	WTTrajListConcat = [[int(i.strip().split()[0]),float(i.strip().split()[1])] for i in WTTrajListConcat]		
	WTTrajListConcat=[i[1] for i in WTTrajListConcat]
	
	#MTHandling
	fileMT = open(fileMT, "r")
	MTTrajListConcat = fileMT.readlines()
	MTTrajListConcat = [[int(i.strip().split()[0]),float(i.strip().split()[1])] for i in MTTrajListConcat]	
	MTTrajListConcat=[i[1] for i in MTTrajListConcat]
	
	# Merge data to find min max value
	merged_list = WTTrajListConcat + MTTrajListConcat
	MinValue=min(merged_list)
	MaxValue=max(merged_list)	
	
	# Plotting
	plt.figure(figsize=(8, 4))
	
	WTweights = np.ones_like(WTTrajListConcat) / len(WTTrajListConcat)
	MTweights = np.ones_like(MTTrajListConcat) / len(MTTrajListConcat)
		
	nWT,binsWT,patchesWT=plt.hist(WTTrajListConcat,bins=100,
				alpha=0.3,color='blue',weights=WTweights,
				label='WT',
				range=[MinValue-(MinValue*0.05), MaxValue+(MaxValue*0.05)])
				
	nMT,binsMT,patchesMT=plt.hist(MTTrajListConcat,bins=100,
				alpha=0.7,color='red',weights=MTweights,
				label='MT',
				range=[MinValue-(MinValue*0.05), MaxValue+(MaxValue*0.05)])	
				
	# Compute the center of each bin for plotting
	bin_centersWT = (binsWT[1:] + binsWT[:-1]) / 2
	bin_centersMT = (binsMT[1:] + binsMT[:-1]) / 2				
	
	# Initial value bar height
	InitialValueHeight=0.0	
	
	if max(nWT) > max(nMT):
		InitialValueHeight=max(nWT)
	else:
		InitialValueHeight=max(nMT)	
	
	InitialValueHeight=InitialValueHeight+ (InitialValueHeight*0.2)
	
	# Initial value bar width
	InitialValueWidth=0.0
	
	widthWT=binsWT[1]-binsWT[0]
	widthMT=binsMT[1]-binsMT[0]		
	
	if widthWT > widthMT:
		InitialValueWidth=widthWT
	else:
		InitialValueWidth=widthMT	

	#print(InitialValueWidth)	
	InitialValueWidth/=2

	# Plotting		
	plt.plot(bin_centersWT, nWT, linewidth=2, color='blue')
	plt.plot(bin_centersMT, nMT, linewidth=2, color='red')
	
	plt.bar([WTInitialValue],[InitialValueHeight],color='green',width=InitialValueWidth,
		label='WT starting volume value')
	plt.bar([MTInitialValue],[InitialValueHeight],color='orange',width=InitialValueWidth,
		label='MT starting volume value')
	plt.tight_layout()	
	plt.xlabel('Volume ($\mathrm{\AA^3}$)',labelpad=10,fontsize=15)
	plt.ylabel('Probability distribution',labelpad=10,fontsize=15)
	plt.legend(loc='upper right')
	plt.savefig(f"{FileOut}.png",dpi=1200,bbox_inches='tight')
	plt.savefig(f"{FileOut}.eps",dpi=1200,bbox_inches='tight')
