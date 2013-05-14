This file describes the experiments performed to validate the estimation of the TRE.  

The main script is masterTRESimulations.m which is actually a function with lots of command line options.  Use the help to see all the command line options.

There are four main experiments run by this script:

1. Test the TRE estimation formula for anisotropic FLE with random rigid fiducial configurations. Published in thesis and IEEE TMI 2008 paper.
	
	Matlab command to recreate data in thesis and paper:
		
		masterTRESimulations('DataPath', '<your choice>', 'RandomFiducial', 1, 'NumSamples', 10^5, 'NumTrials', 1000, 'NumOrientations', 100 );
	
	The other default values as described in help masterTRESimulations are maintained.
		
2. Test the TRE estimation formula for anisotropic FLE with a particular probe rigid body configuration. Published in thesis and IEEE TMI 2008 paper.
	
	Matlab command to recreate data in thesis and paper:
		
		masterTRESimulations('DataPath', '<your choice>', 'RandomFiducial', 0,'RandomProbe', 1, 'SigmaIso', diag([1/27, 1/27, 1/27]), 'SigmaAniso', diag([1/99, 1/99, 1/11]), 'NumSamples', 10^5, 'NumOrientations', 1000, 'WestProbeDesign', {'d', 'e'}, 'WestProbe_rho', [85, 170]);
	
	The other default values as described in help masterTRESimulations are maintained.
		
3. Test the TRE estimation formula for anisotropic FLE with a particular probe and reference tool configuration.  Published in MICCAI 2007 paper.
	
	Matlab command to recreate data in thesis and paper:
		
		masterTRESimulations('DataPath', '<your choice>', 'RandomFiducial', 0,'RandomProbeRef', 1, 'SigmaIso', diag([1/27, 1/27, 1/27]), 'SigmaAniso', diag([1/99, 1/99, 1/11]), 'NumSamples', 10^5, 'NumOrientations', 1000, 'WestRef_r', [32, 64]); 
		
4. Test the TRE estimation formula for non-homogenous anisotropic FLE with random rigid fiducial configurations.
	
	Matlab command to recreate data in thesis:
		
		1. FLE Model 1 - radial model, see getFLEMatrix.m or pg. 112 of Wiles Thesis.
			masterTRESimulations('DataPath', '<your choice>', 'RandomFiducial', 0,'RandomFiducialNH', 1, 'NumBodies', 100, 'NumTrials', 1000, 'NumOrientations', 10, 'NH_FLEModelType', 1, 'NH_FLEWeight', [1.0, 1.5, 3.0])
		
		2. FLE Model 2 - random model, see getFLEMatrix or pg. 112 of Wiles Thesis.
			masterTRESimulations('DataPath', '<your choice>', 'RandomFiducial', 0,'RandomFiducialNH', 1, 'NumBodies', 100, 'NumTrials', 1000, 'NumOrientations', 10, 'NH_FLEModelType', 2);
			
		3. FLE Model 3 - random model and rotated FLE, see getFLEMatrix or pg. 112 of Wiles Thesis.
			masterTRESimulations('DataPath', '<your choice>', 'RandomFiducial', 0,'RandomFiducialNH', 1, 'NumBodies', 100, 'NumTrials', 1000, 'NumOrientations', 10, 'NH_FLEModelType', 3);
