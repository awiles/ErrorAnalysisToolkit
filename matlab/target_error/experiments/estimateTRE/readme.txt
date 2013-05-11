This file describes the experiments performed to validate the estimation of the TRE.  

The main script is masterTRESimulations.m which is actually a function with lots of command line options.  Use the help to see all the command line options.

There are four main experiments run by this script:

1. Test the TRE estimation formula for anisotropic FLE with random rigid fiducial configurations. Published in thesis and IEEE TMI 2008 paper.
	Matlab command to recreate data in thesis and paper:
		masterTRESimulations('DataPath', '<your choice>', ...
		
2. Test the TRE estimation formula for anisotropic FLE with a particular probe rigid body configuration. Published in thesis and IEEE TMI 2008 paper.
	Matlab command to recreate data in thesis and paper:
		masterTRESimulations('DataPath', '<your choice>', ...
		
3. Test the TRE estimation formula for anisotropic FLE with a particular probe and reference tool configuration.  Published in MICCAI 2007 paper.
	Matlab command to recreate data in thesis and paper:
		masterTRESimulations('DataPath', '<your choice>', ...
		
4. Test the TRE estimation formula for non-homogenous anisotropic FLE with random rigid fiducial configurations.
	Matlab command to recreate data in thesis:
		masterTRESimulations('DataPath', '<your choice>', ...
