function masterTRESimulations(varargin)
% masterTRESimulations.m
%**************************************************************************
% Master Simulation Script to test TRE formulations.
%
%**************************************************************************
%   This script will perform Monte Carlo simulations for different
%   simulation test cases and compare the Monte Carlo results to the
%   theoretical results provided by the various formulations presented in
%   my thesis.
%   
%   There are several optional variables that can be passed in.  If nothing
%   is passed in then the basic simulation for Random Rigid bodies with the
%   default parameters below.
%   
%   Optional Parameters:
%       'DataPath'      High-level path to write the data and results.
%                           - default: userpath() <-- OS dependent.
%       'RandomFiducial' Run the random fiducial experiment.
%                           - default: 1 (On)
%                           - IEEE TMI March 2008 - Random Configurations.
%       'RandomProbe'   Run the random probe experiment.
%                           - default: 0 (Off)
%                           - IEEE TMI March 2008 - West Tools.
%       'RandomProbeRef' Run the random probe experiment w/ reference tool.
%                           - default: 0 (Off)
%                           - MICCAI 2007.
%       'RandomFiducialNH' Run the random fiducial non-homogenous FLE experiment.
%                           - default: 0 (Off)
%                           - Published in thesis only.
%       'SigmaIso'      Pass in a custom isotropic FLE 3 x 3 matrix.
%                           - default: (1/3)*eye(3), FLE_RMS = 1mm
%       'SigmaAniso'    Pass in a custom anistorpic FLE 3 x 3 matrix.
%                           - default: diag([1/11, 1/11, 9/11]), FLE_RMS = 1mm
%       'NumSamples'    Number of random samples in each MC simulation.
%                           - default: 10^4
%       'NumBodies'     Number of random rigid bodies to be tested.
%                           - default: 10
%                           - only used in the random rigid bodies.
%       'NumTrials'     Number of trials for each rigid body configuration.
%                           - default: 1
%       'NumOrientations' Number of random orientations tested.
%                           - default: 10
%       'NumPositions'  Number of positions tested - non-homogenous only.
%                           - default: 1
%       'NumMarkers'    Number of markers in random fiducial tests.
%                           - default: [3 4 10 20 50]
%       'MarkerRange'   Range over which the random fiducials can be selected.
%                           - default: 200
%       'TargetRange'   Range over which the random targets can be selected.
%                           - default: 400
%       'WestProbeDesign' Probe configuration from West paper, circa 2004.
%                           - options: 'd' or 'e' <-- can be a cell array.
%                           - default: 'e'
%       'WestProbe_A'   Pointer Probe dimension, parameter A - height.
%                           - default: 71
%       'WestProbe_B'   Pointer Probe dimension, parameter B - width.
%                           - default: 54
%       'WestProbe_rho' Pointer Probe dimension, parameter rho - tool tip offset.
%                           - default: 85 <-- this can be a vector.
%       'WestRef_r'     Reference tool dimension, parameter r - square size.
%                           - default: 32 <-- this can be a vector.
%       'WestRef_d'     Reference tool distance to Pointer Probe tool tip.
%                           - default: [100, 200, 300, 400]
%       'NH_FLEModelType' FLE model type for non-homogenous FLE.
%                           - options: TODO
%                           - default: 1 - spherical model
%       'NH_FLEOrigin'  FLE model origin.
%                           - default: [0 0 0]
%       'NH_FLERate'    FLE model rate of increase from origin in mmRMS/mm.
%                           - default: 0.002
%       'NH_FLEBaseline' FLE RMS at origin.
%                           - default: 0.1
%       'NH_FLEWeight'  FLE relative weights for each axis. Still assuming indepenent axis values.
%                           - default: [1 1 3]
%**************************************************************************

%% defaults.

% high level path:
datapath = userpath();

% list of experiments to process.
bRandFiducial = 1;  % IEEE TMI March 2008 - Random Configurations.
bRandProbe = 0;     % IEEE TMI March 2008 - West Tools.
bRandProbeRef = 0;  % MICCAI 2007.
bRandFiducialNH = 0;% Non-Homogenous simulations.

% FLE_RMS = 1 mm
SigmaIso = (1/3)*eye(3);
SigmaAniso = diag([1/11, 1/11, 9/11]);

% test sample parameters.
nSamples = 10^4;
nBodies = 10;
nTrials = 1;
nOrientations = 10;
nPositions = 1;

% parameters for the random fiducial tests.
nMarkers = [3 4 10 20 50];
mrkRange = 200;
targetRange = 400;

% probe design testing.
% West tool design parameters.
design = {'e'};
A = 71;
B = 54;
rho = 85;
% reference parameters.
r = 32;
d = [100, 200, 300, 400];

% non-homogenous FLE parameters.
fleparms.modeltype = 1;     %spherical model.
fleparms.origin = [0 0 0];  %origin of sphere.
fleparms.rate = 0.002;      %rate of increase from origin in mmRMS/mm
fleparms.baseline = 0.1;    %RMS at origin.
fleparms.weight = [1 1 3];  %relative weights for each axis. Still assuming indepenent axis values.

%% get the variable arguments if any:
if( nargin > 0 )
    nVarArgs = length(varargin);
    i = 1;
    while( i <= nVarArgs )
        if( strcmp(varargin{i}, 'DataPath' ) )
            i = i+1;
            datapath = varargin{i};
            % check that the path exists.
            if( ~isdir(datapath) )
                error( 'Invalid path given: %s', datapath);
            end
        elseif( strcmp(varargin{i}, 'RandomFiducial') )
            i = i+1;
            bRandFiducial = varargin{i};
        elseif( strcmp(varargin{i}, 'RandomProbe') )
            i = i+1;
            bRandProbe = varargin{i};
        elseif( strcmp(varargin{i}, 'RandomProbeRef') )
            i = i+1;
            bRandProbeRef = varargin{i};
        elseif( strcmp(varargin{i}, 'RandomFiducialNH') )
            i = i+1;
            bRandFiducialNH = varargin{i};
        elseif( strcmp(varargin{i}, 'SigmaIso') )
            i = i+1;
            SigmaIso = varargin{i};
            if( size(SigmaIso ~= [3 3]) )
                error('Invalid matrix given for SigmaIso, please check your inputs.');
            end
        elseif( strcmp(varargin{i}, 'SigmaAniso') )
            i = i+1;
            SigmaAniso = varargin{i};
            if( size(SigmaIso ~= [3 3]) )
                error('Invalid matrix given for SigmaAniso, please check your inputs.');
            end
        elseif( strcmp(varargin{i}, 'NumSamples') )
            i = i+1;
            nSamples = varargin{i};
        elseif( strcmp(varargin{i}, 'NumBodies') )
            i = i+1;
            nBodies = varargin{i};
        elseif( strcmp(varargin{i}, 'NumTrials') )
            i = i+1;
            nTrials = varargin{i};
        elseif( strcmp(varargin{i}, 'NumOrientations') )
            i = i+1;
            nOrientations = varargin{i};
        elseif( strcmp(varargin{i}, 'NumPositions') )
            i = i+1;
            nPositions= varargin{i};
        elseif( strcmp(varargin{i}, 'NumMarkers') )
            i = i+1;
            nMarkers = varargin{i};
        elseif( strcmp(varargin{i}, 'MarkerRange') )
            i = i+1;
            mrkRange = varargin{i};
        elseif( strcmp(varargin{i}, 'TargetRange') )
            i = i+1;
            targetRange = varargin{i};
        elseif( strcmp(varargin{i}, 'WestProbeDesign') )
            i = i+1;
            if( iscell(varargin{i}) )
                design = varargin{i};
            elseif( ischar(varargin{i}) )
                design = varargin(i);
            else
                error('Invalid value given for WestProbeDesign.');
            end
        elseif( strcmp(varargin{i}, 'WestProbe_A') )
            i = i+1;
            A = varargin{i};
        elseif( strcmp(varargin{i}, 'WestProbe_B') )
            i = i+1;
            B = varargin{i};
        elseif( strcmp(varargin{i}, 'WestProbe_rho') )
            i = i+1;
            rho = varargin{i};
        elseif( strcmp(varargin{i}, 'WestRef_r') )
            i = i+1;
            r = varargin{i};
        elseif( strcmp(varargin{i}, 'WestRef_d') )
            i = i+1;
            d = varargin{i};
        elseif( strcmp(varargin{i}, 'NH_FLEModelType') )
            i = i+1;
            fleparms.modeltype = varargin{i};
        elseif( strcmp(varargin{i}, 'NH_FLEOrigin') )
            i = i+1;
            fleparms.origin = varargin{i};
        elseif( strcmp(varargin{i}, 'NH_FLERate') )
            i = i+1;
            fleparms.rate = varargin{i};
        elseif( strcmp(varargin{i}, 'NH_FLEBaseline') )
            i = i+1;
            fleparms.baseline = varargin{i};
        elseif( strcmp(varargin{i}, 'NH_FLEWeight') )
            i = i+1;
            fleparms.weight = varargin{i};
        else
            error('Unknown parameter: %s', varargin{i})
        end
        i=i+1;
    end
end


%% Move into the correct high level path, if it exists.
if( isdir(datapath) )
    cd(datapath)
    if( ~isdir('TRESimulations') )
        mkdir('TRESimulations')
    end
    cd('TRESimulations');
else
    error( 'Invalid path given: %s', datapath);
end

%% Set up the FLE options.
% TODO: make iso optional.
Sigma = {SigmaIso SigmaAniso};

%% Monte Carlo simulation for random fiducials.
% Published in IEEE TMI March 2008 - Random Configurations.
if(bRandFiducial)
    if( ~isdir('RandomFiducial') )
        mkdir('RandomFiducial');
    end
    cd('RandomFiducial');
    
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(Sigma)
        for j = 1:length(nMarkers)
            nCount = nCount + 1;
            if(nCount > nDone)
                runRandFiducialTargetExperiment(nSamples, nBodies,...
                    nTrials, nOrientations, nPositions,...
                    nMarkers(j), Sigma{i}, mrkRange, targetRange);
            end
        end
    end
    cd ..
end

%% Monte Carlo simulation for probe in multiple orientations.
% Published in IEEE TMI March 2008 - West Tools.
if(bRandProbe)
    if( ~isdir('RandomProbe') )
        mkdir('RandomProbe');
    end
    cd('RandomProbe');
    
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(Sigma)
        for j = 1:length(design)
            for k = 1:length(rho)
                nCount = nCount + 1;
                if(nCount > nDone)
                    runRandProbeXfrmExperiment(nSamples,...
                        nOrientations, nPositions, design{j}, A, B, rho(k), Sigma{i});
                end
            end
        end
    end
    cd ..
end


%% Monte Carlo simulation for probe being tracked relative to a reference.
% Published in MICCAI 2007.
if (bRandProbeRef)
    if( ~isdir('RandomProbeRef') )
        mkdir('RandomProbeRef');
    end
    cd('RandomProbeRef');
    
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(Sigma)
        for j = 1:length(design)
            for k = 1:length(rho)
                for m = 1:length(r)
                    for n = 1:length(d)
                        nCount = nCount + 1;
                        if(nCount > nDone)
                            runRandProbeRefXfrmExperiment(nSamples,...
                                nOrientations, nPositions, design{j}, ...
                                A, B, rho(k), r(m), d(n), Sigma{i});
                        end
                    end
                end
            end
        end
    end
    cd ..
end

%% Monte Carlo simulation for random fiducials and non-homogenous FLE .
% Not published (thesis only) - Random Configurations.
if(bRandFiducialNH)
    if( ~isdir('RandomFiducialNH') )
        mkdir('RandomFiducialNH');
    end
    cd('RandomFiducialNH');
    
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(fleparms)
        for j = 1:length(nMarkers)
            nCount = nCount + 1;
            if(nCount > nDone)
                runRandFiducialTargetNonHomogenousExperiment...
                    (nSamples, nBodies, nTrials, nOrientations, nPositions,...
                    nMarkers(j), fleparms, mrkRange, targetRange)
            end
        end
    end
    cd ..
end

% return to main data path.
cd(datapath)