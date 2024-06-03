The repository contains the results of molecular dynamics simulations conducted in the publication "Theoretical Investigations of a point mutation on the H5N1 Hemagglutinin-Receptor Complex."

	. Directory "PDBStructure_Topology" contains the WT and MT D94N H5 Hemaggluinin structures, and SAα2,3Gal and SAα2,6Gal glycans structure and topologies in Glycam and Charmm formats.
	
	. Directory "MDSimulationScripts" contains computational-workflow to conduct (1) complex of HA-glycans, (2) free glycans multi-simulations, and (3) Binding free energy change of HA to the glycans upon the D94N mutation.
		
	. Directory "MultiSimulationsAnalysis" includes detailed analyses for the simulations of HA-glycans complexes and free glycans. 
		. Sub-directory "PlottingScripts" includes the python scripts that were used to generate the Figures in the publication.
		. "AnalysisScript_WT_MTHA_SA23.sh": automation script to conduct the analyses on the complex of WT/MT HA with glycan SA23
		. "AnalysisScript_WT_MTHA_SA26.sh": automation script to conduct the analyses on the complex of WT/MT HA with glycan SA26
		
		For complex HA-glycans simulations
			. Sub-directory "ProteinReceptorDistance" contains the data to generate Figure 3.
			. Sub-directory "DetectingHbondTargetProteinRes_System" contains the data to generate Figure 4, S4.
			. Sub-directory "DetectingHbondReceptor_Protein" contains the data to generate Figure 5.
			. Sub-directories "ReceptorRMSD" and "Angle" contain the data to generate Figure 8.
			. Sub-directories "ProteinReceptorDistance" and "Energy" contain the data to generate Figure S3.
			. Sub-directories "RMSF" contain the data to generate Figure S5.
			
		For free glycans simulations
			
			. Sub-directories "ReceptorRMSD" and "Angle" contain the data to generate Figure 7.
