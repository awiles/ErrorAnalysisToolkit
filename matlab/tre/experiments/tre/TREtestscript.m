% TREtestscript

cd 'e:\data\tretest\'

%% FLE Models.
% 1. FLE_RMS = 1 mm
SigmaIso = (1/3)*eye(3);
SigmaAniso = diag([1/11, 1/11, 9/11]);

% 2. FLE_RMS = 1/3mm
%SigmaIso = (((1/3)^2)/3)*eye(3);
%SigmaAniso = diag([1/99, 1/99, 1/11]);

%% test sample parameters.
nSamples = 10^4;
nBodies = 10;
nTrials = 1;
nOrientations = 10;
nPositions = 1;

%% probe design testing.
% West tool design parameters.
design = {'e'};
A = 71;
B = 54;
rho = [85];
r = [32];
d = [100, 200, 300, 400];

Sigma = {SigmaAniso};

% parameters for the random fiducial tests.
nMarkers = [3 4 10 20 50];
%nMarkers = 4;
mrkRange = 200;
targetRange = 400;

% non-homogenous FLE parameters.
fleparms.modeltype = 1;     %spherical model.
fleparms.origin = [0 0 0];  %origin of sphere.
fleparms.rate = 0.002;      %rate of increase from origin in mmRMS/mm
fleparms.baseline = 0.1;    %RMS at origin.
fleparms.weight = [1 1 3];  %relative weights for each axis. Still assuming indepenent axis values.


%% list of experiments to process.
bRandFiducial = 1;  % IEEE TMI March 2008 - Random Configurations.
bRandProbe = 0;     % IEEE TMI March 2008 - West Tools.
bRandProbeRef = 0;  % MICCAI 2007.
bRandFiducialNH = 0;% Non-Homogenous simulations.

%% Monte Carlo simulation for random fiducials and non-homogenous FLE .
%% Not published - Random Configurations.
if(bRandFiducialNH)
    cd '.\non-homog-testing'
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(fleparms)
        for j = 1:length(nMarkers)
            nCount = nCount + 1;
            if(nCount > nDone)
                [probe,rmsDiff] = runRandFiducialTargetNonHomogenousExperiment...
                    (nSamples, nBodies, nTrials, nOrientations, nPositions,...
                    nMarkers(j), fleparms, mrkRange, targetRange)
            end
        end
    end
end

%% Monte Carlo simulation for random fiducials.
%% Published in IEEE TMI March 2008 - Random Configurations.
if(bRandFiducial)
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(Sigma)
        for j = 1:length(nMarkers)
            nCount = nCount + 1;
            if(nCount > nDone)
                [probe,rmsDiff] = runRandFiducialTargetExperiment(nSamples, nBodies,...
                    nTrials, nOrientations, nPositions,...
                    nMarkers(j), Sigma{i}, mrkRange, targetRange)
            end
        end
    end
end

%% Monte Carlo simulation for probe in multiple orientations.
%% Published in IEEE TMI March 2008 - West Tools.
if(bRandProbe)
    % quick variables to start midway through.
    nCount = 0;
    nDone = 0;
    for i = 1:length(Sigma)
        for j = 1:length(design)
            for k = 1:length(rho)
                nCount = nCount + 1;
                if(nCount > nDone)
                    [probe,rmsDiff]=runRandProbeXfrmExperiment(nSamples,...
                        nOrientations, nPositions, design{j}, A, B, rho(k), Sigma{i});
                end
            end
        end
    end
end


%% Monte Carlo simulation for probe being tracked relative to a reference.
%% Published in MICCAI 2007.
if (bRandProbeRef)
    quick variables to start midway through.
    nCount = 0;
    nDone = 6;
    for i = 1:length(Sigma)
        for j = 1:length(design)
            for k = 1:length(rho)
                for m = 1:length(r)
                    for n = 1:length(d)
                        nCount = nCount + 1;
                        if(nCount > nDone)
                            [probe,rmsDiff]=runRandProbeRefXfrmExperiment(nSamples,...
                                nOrientations, nPositions, design{j}, A, B, rho(k), r(m), d(n), Sigma{i});
                        end
                    end
                end
            end
        end
    end
end