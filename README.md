# The repository contains molecular dynamics simulations data presented in the publication "Theoretical Investigations of a point mutation affecting H5 Hemagglutinin’s receptor binding preference."

## 1. Research summary:
 In this research, we aim to elucidate the potential impact of a single amino acid substitution, D94N on the surface of the influenza virus hemagglutinin (HA) protein. HA is a glycoprotein that is essential for viral entry into host cells through recognizing and interacting with the host Sialic Acid (SA) containing receptor. The HA from influenza virus strains circulating in human and avian species exhibit a preference for interacting with the receptor containing SA linked to galactose via an α-2,6 bond (SAα2,6Gal) and via an α-2,3 bond (SAα2,3Gal), respectively. Thus, the SAα2,6Gal and SAα2,3Gal are commonly referred to as the “human receptor” and “avian receptor”, respectively. 
 
 Even though the wild-type amino acid D94 does not directly interact with the glycan receptor, Previous experimental research demonstrated that the D94N mutation induces a switch in HA receptor binding preference from the avian-type SAα2,3Gal to the human-type SAα2,6Gal. Our computational research, employing molecular dynamics simulations and binding free energy computations, corroborated this change in receptor preference. Specifically, the mutation was observed to cause reduced stability for the HA-avian-type SAα2,3Gal complex and enhanced stability for the HA-human-type SAα2,6Gal. Our subsequent analysis clearly illustrated that the mutation's detrimental effect on the HA-SAα2,3Gal interactions is due to increased flexibility of the 130-Loop, which is crucial for glycan recognition by HA. Conversely, the resilience in the interactions between HA and SAα2,6Gal is attributed to the inherent flexibility of the glycan.
 
 Our research presents a simple but effective strategy to enhance the sampling of the HA-glycan complexes (or protein-small molecule in general). The Python-based in-house design analyses, including binding probability analysis, HA receptor binding site volume computation, glycan topology classification, and particularly ‘hydrogen bond formation propensity’ analysis, has proven effective in deciphering the underlying molecular mechanisms and sequential effects of the amino acid substitution. To conclude, our simulation strategy and analyses can be valuable for predicting and deciphering phenomena that are difficult to capture through wet-lab experiments or to probe the interactions between small molecules/glycans and their biological targets. 
 
 Finally, to read the full publication, please access the following link https://www.sciencedirect.com/science/article/pii/S1476927124001774

## 2.Repository structure

Directory "PDBStructure_Topology" contains the WT and MT D94N H5 Hemaggluinin structures, and SAα2,3Gal and SAα2,6Gal glycans structure and topologies in Glycam and Charmm formats.

Directory "MDSimulationScripts" contains computational-workflow to conduct (1) complex of HA-glycans, (2) free glycans multi-simulations, and (3) Binding free energy change of HA to the glycans upon the D94N mutation.
		
Directory "MultiSimulationsAnalysis" includes detailed analyses for the simulations of HA-glycans complexes and free glycans. 
	
		Sub-directory "PlottingScripts" includes the python scripts that were used to generate the Figures in the publication.
 		"AnalysisScript_WT_MTHA_SA23.sh": automation script to conduct the analyses on the complex of WT/MT HA with glycan SA23
   		"AnalysisScript_WT_MTHA_SA26.sh": automation script to conduct the analyses on the complex of WT/MT HA with glycan SA26
     
     		For each sub-directory ("WTHA_SA23", "MTHA_SA23", "WTHA_SA26", "MTHA_SA26") containing data obtained from complex of HA-glycans simulations, the sub-directory is organized as follows: 
       			. Directory "ProteinReceptorDistance" contains the data to generate Figure 3 of the publication.
			. Directory "DetectingHbondTargetProteinRes_System" contains the data to generate Figure 4, S4 of the publication.
   			. Directory "DetectingHbondReceptor_Protein" contains the data to generate Figure 5 of the publication.
      			. Directories "ReceptorRMSD" and "Angle" contain the data to generate Figure 8 of the publication.
	 		. Directories "ProteinReceptorDistance" and "Energy" contain the data to generate Figure S3 of the publication.
    			. Directories "RMSF" contain the data to generate Figure S5 of the publication.
       			. Directories "RBS_volume_calculation" contain the data to generate Figure 9 of the publication.
	  
   		The two directories "FreeSA23" and "FreeSA26" contain data obtained from free glycans in solution simulations
     			. Directories "ReceptorRMSD" and "Angle" contain the data to generate Figure 7 of the publication.


To re-produce any analysis described in the publication for user own research system, please consulting either the bash script "AnalysisScript_WT_MTHA_SA23.sh" or "AnalysisScript_WT_MTHA_SA26.sh". To generate plots or graphs using the raw data achieved from our designed analysis workflow, considering utilizing the python scripts in sub-directory directory 'MultiSimulationsAnalysis/PlottingScripts'
    
## 3. Additional information

To create environment for MDAnalysis: conda env create -f env.yml

Directory 'MultiSimulationsAnalysis/PlottingScripts' contains the python scripts to generate figures presented in the publication. Invoking 'python3 file_name --help' to access the description of the program and instruction, how to use the program.

The script in MDAnalysis directory requires software and python libraries specified in env.yml script.
