function [probe,rmsDiff] = runRandFiducialTargetNonHomogenousExperiment(nSamples,...
    nBodies, nTrials, nOrientations, nPositions,nMarkers, fleparms, mrkRange, targetRange, varargin)
% [probe,rmsDiff] = runRandFiducialTargetNonHomogenousExperiment(nSamples, nTrials,...
%                       nOrientations, nPositions, nMarkers, Sigma, mrkRange, targetRange)
%
%       nSamples - number of random samples in each Monte Carlo simulation,
%                   i.e. N = 10,000 or 10^6 or something similar.
%       nBodies  - number of random rigid body configurations that are to be tested.
%       nTrials  - number of Monte Carlo simulations for for each rigid
%                   body configuration.
%       nOrientations - number of random orientations to test for each
%                   trial.
%       nPositions - number of random positions to test (N/A if Sigma is
%                       constant throughout the volume, therefore set to
%                       one.  Kept here as a place holder for future work.)
%       nMarkers - number of fiducial markers randomly chosen for each
%                       body trial.
%       fleparms    - parameters for the non-homogenous FLE model.
%       mrkRange - distance over which the markers will appear.  The origin
%                       is assumed to bisect this range.
%       targetRange - distance over which the target will appear.  The
%                       origin is assumed to bisect this range.

% defaults for optional arguments.
smtp_server = '';
email_address = '';
bEmail = 0;

if( nargin > 9 )
    nVarArgs = length(varargin);
    i = 1;
    while( i <= nVarArgs )
        if( strcmp(varargin{i}, 'smtp') )
            i=i+1;
            smtp_server = varargin{i};
        elseif (strcomp(varargin{i}, 'email'))
            i=i+1;
            email_address = varargin{i};
        elseif( strcmp(varargin{i}, 'Verbose'))
            verbose = 1;
        else
            error('Unknown paramter: %s', varargin{i});
        end
    end
end

if( ~isempty(smtp_server) && ~isempty(email_address))
    email_setup('smtp.me.com', 'awiles@me.com');
    bEmail = 1;
end

nCount=0;
nTotalCount = nBodies*nTrials*nOrientations*nPositions;
passCount1 = 0;
passCount2 = 0;
passCount3 = 0;
passCount4 = 0;
rmsDiff = zeros(nTotalCount,4);
rmsDiffAvgFLE = rmsDiff;

% set up the data directory to store the results.
starttime = clock;
datadir = sprintf('%02dMarkers-%03dBodies-%03dOrientations-%03dPositions-%02dFLEModel',...
    nMarkers, nBodies, nOrientations, nPositions, fleparms.modeltype);
if( ~isdir(datadir) )
    mkdir(datadir);
end
cd(datadir);
datetime   = sprintf('%4d%02d%02d-%02d%02d', starttime(1:5));
if( ~isdir(datetime) )
    mkdir(datetime);
end
cd(datetime);

% write the experiment parameters to the readme.txt file.
frm = fopen('readme.txt', 'wt');
fprintf(frm, 'Non-Homogenous FLE Model testing\n');
fprintf(frm, 'Date processed: %s\n', datetime); parm.name = datetime;
fprintf(frm, 'nSamples = %d\n', nSamples); parm.nSamples = nSamples;
fprintf(frm, 'nBodies = %d\n', nBodies); parm.nBodies = nBodies;
fprintf(frm, 'nTrials = %d\n', nTrials); parm.nTrials = nTrials;
fprintf(frm, 'nOrientations = %d\n', nOrientations); parm.nOrientations = nOrientations;
fprintf(frm, 'nPositions = %d\n', nPositions); parm.nPositions = nPositions;
fprintf(frm, 'nMarkers = %d\n', nMarkers); parm.nMarkers = nMarkers;
fprintf(frm, 'mrkRange = %d\n', mrkRange); parm.mrkRange = mrkRange;
fprintf(frm, 'targetRange = %d\n\n', targetRange); parm.targetRange = targetRange;
%TODO: update here:
fprintf(frm, 'FLE Parms \n');
fprintf(frm, 'modeltype = %d\n', fleparms.modeltype); parm.fleparms = fleparms;
fprintf(frm, 'origin = [%f, %f, %f]\n', fleparms.origin);
fprintf(frm, 'rate = %f\n', fleparms.rate);
fprintf(frm, 'baseline = %f\n', fleparms.baseline);
fprintf(frm, 'weight = [%f, %f, %f]\n\n', fleparms.weight);

fprintf(frm, 'Data Summary\n------------\n\n');

% set up the data file.
fdata = fopen('data.csv', 'wt');

% perform the simulation.
for i=1:nBodies
    % build the marker set.
    probe = generateRandFiducialTargets(nMarkers, mrkRange, targetRange);
    probe.bHomogenous = 0;
    % DEBUG:
    %     A = 71;
    %     B = 54;
    %     rho = 85;
    %     probe.Rigid.mrk = [ 0 B/2 0;
    %         -A/2 0 0;
    %         0 -B/2 0;
    %         A/2 0 0 ];
    %     probe.Rigid.normals = repmat([ 0 0 1], 4, 1);
    %     probe.Rigid.tip = [ -rho 0 0 ];
    % ENDDEBUG



    for j=1:nTrials
        for k=1:nOrientations
            %******************Get the Random Orientation *******%
            [probeXfrm.R, q] = getRandOrientation();
            for m=1:nPositions
                %*********************Get the Random Position *********%
                %assuming no translation for now.
                probeXfrm.pos = [0 0 0];
                %****************Transform into the Random Xfrm*******%
                %move probe coordinate frame to tip.
                %temp.Rigid.mrk = probe.Rigid.mrk - repmat(probe.Rigid.tip, size(probe.Rigid.mrk,1), 1);
                %temp.Rigid.tip = [0 0 0];
                %transform probe rigid body into test space.
                probe.Actual.mrk = (probeXfrm.R * probe.Rigid.mrk')'...
                    + repmat(probeXfrm.pos, size(probe.Rigid.mrk,1), 1);
                probe.Actual.tip = (probeXfrm.R * probe.Rigid.tip')' + probeXfrm.pos;

                %**************Compute the theoretical statistics.*******%
                % set-up the non-homogenous model parameters.
                sigmaNH = zeros(3,3,nMarkers);
                for n = 1:nMarkers
                    sigmaNH(:,:,n) = getFLEMatrix(fleparms, probe.Actual.mrk(n,:));
                end
                %sigmaNH;
                sigmaAvg = mean(sigmaNH,3);
                probe.FLE.sigmaNH = sigmaNH;
                probe.FLE.sigmaAvg = sigmaAvg;
                probe.Actual.stats.mu = zeros(1,3);
                probe.AvgFLE.stats.mu = zeros(1,3);
                %use the non-homogenous model.
                %fprintf('Get NH TRE...');
                [probe.Actual.stats.RMS, probe.Actual.stats.Sigma, probe.Actual.stats.SigmaPA]...
                    = calcTRE(sigmaNH, [probe.Actual.mrk; probe.Actual.tip]);
                %[probe.Actual.mrk; probe.Actual.tip]
                %fprintf('NH   RMS: %f\n', probe.Actual.stats.RMS);   
                %average the FLE matrices and use the homogenous model.
                %fprintf('Get Avg TRE...\n');
                [probe.AvgFLE.stats.RMS, probe.AvgFLE.stats.Sigma, probe.AvgFLE.stats.SigmaPA]...
                    = calcTRE(sigmaAvg, [probe.Actual.mrk; probe.Actual.tip]);
                %fprintf('Avg  RMS: %f\n', probe.AvgFLE.stats.RMS);
            
                %estimate using Fitzpatrick's.
                probe.Fitz.stats.RMS = calcTREFitz(sqrt(trace(sigmaAvg)), [probe.Actual.mrk; probe.Actual.tip]);
                %fprintf('Fitz RMS: %f\n', probe.Fitz.stats.RMS);
                %fprintf('NH   Sigma: \n');
                %sigmaNH
                %fprintf('Avg  Sigma: \n');
                %sigmaAvg
                
                %*******************Monte Carlo trial.**********%
                [probe.Meas.error, probe.Meas.tip] = simTRE(sigmaNH, nSamples, probe.Rigid, probe.Actual );
                probe = computeStats(probe);
                %compare the simulation to the theoretical stats.
                testResult = compareTREStats(probe);
                passCount1 = passCount1 + testResult.probe.covariance;
                passCount2 = passCount2 + testResult.probe.meanandcov;
                passCount3 = passCount3 + testResult.probeAvgFLE.covariance;
                passCount4 = passCount4 + testResult.probeAvgFLE.meanandcov;
                nCount = nCount + 1;
                rmsDiff(nCount,1) = testResult.probe.RMS.Meas;
                rmsDiff(nCount,2) = testResult.probe.RMS.Theory;
                rmsDiff(nCount,3) = testResult.probe.RMS.Diff;
                rmsDiff(nCount,4) = testResult.probe.RMS.PercentDiff;
                rmsDiffAvgFLE(nCount,1) = testResult.probeAvgFLE.RMS.Meas;
                rmsDiffAvgFLE(nCount,2) = testResult.probeAvgFLE.RMS.Theory;
                rmsDiffAvgFLE(nCount,3) = testResult.probeAvgFLE.RMS.Diff;
                rmsDiffAvgFLE(nCount,4) = testResult.probeAvgFLE.RMS.PercentDiff;

                fprintf('Completed %3.1f%%, nMarkers = %d, RMS - NH: %3.4f (%3.2f%%), Avg: %3.4f (%3.2f%%), MC: %3.4f, Fitz: %3.4f  \n',...
                    (100*nCount/nTotalCount), nMarkers, probe.Actual.stats.RMS, testResult.probe.RMS.PercentDiff, ...
                    probe.AvgFLE.stats.RMS, testResult.probeAvgFLE.RMS.PercentDiff, ...
                    testResult.probe.RMS.Meas, probe.Fitz.stats.RMS);

                % save the data.
                filename = sprintf('data%06d', nCount);
                save(filename, 'probe', 'testResult');
                % write out the summary data.
                fprintf(fdata,'%d, %d, %d, %3.6f, %3.6f, %3.6f, %3.2f, %d, %d, %3.6f, %3.6f, %3.6f, %3.2f\n',...
                    nCount, ...
                    testResult.probe.covariance, testResult.probe.meanandcov,...
                    testResult.probe.RMS.Meas, testResult.probe.RMS.Theory,...
                    testResult.probe.RMS.Diff, testResult.probe.RMS.PercentDiff, ...
                    testResult.probeAvgFLE.covariance, testResult.probeAvgFLE.meanandcov,...
                    testResult.probeAvgFLE.RMS.Meas, testResult.probeAvgFLE.RMS.Theory,...
                    testResult.probeAvgFLE.RMS.Diff, testResult.probeAvgFLE.RMS.PercentDiff);
            end
        end
    end
end

%compute the percentage of successful Wishart tests.
passPercent1 = 100*passCount1/nCount;
passPercent2 = 100*passCount2/nCount;
passPercent3 = 100*passCount3/nCount;
passPercent4 = 100*passCount4/nCount;
% compute the stats of the RMS difference.
rmsDiffStats.mu = mean(rmsDiff(:,3));
rmsDiffStats.std = std(rmsDiff(:,3));
rmsDiffStats.max = max(rmsDiff(:,3));
rmsDiffStats.min = min(rmsDiff(:,3));
rmsDiffAvgFLEStats.mu = mean(rmsDiffAvgFLE(:,3));
rmsDiffAvgFLEStats.std = std(rmsDiffAvgFLE(:,3));
rmsDiffAvgFLEStats.max = max(rmsDiffAvgFLE(:,3));
rmsDiffAvgFLEStats.min = min(rmsDiffAvgFLE(:,3));
% compute the stats of the RMS percent difference.
rmsPercentDiffStats.mu = mean(rmsDiff(:,4));
rmsPercentDiffStats.std = std(rmsDiff(:,4));
rmsPercentDiffStats.max = max(rmsDiff(:,4));
rmsPercentDiffStats.min = min(rmsDiff(:,4));
rmsPercentDiffAvgFLEStats.mu = mean(rmsDiffAvgFLE(:,4));
rmsPercentDiffAvgFLEStats.std = std(rmsDiffAvgFLE(:,4));
rmsPercentDiffAvgFLEStats.max = max(rmsDiffAvgFLE(:,4));
rmsPercentDiffAvgFLEStats.min = min(rmsDiffAvgFLE(:,4));

% write summary to readme file.
fprintf(frm, 'Non-Homogenous Model: Wishart Test Results\n-------------------------------\n');
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', passPercent1);
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', passPercent2);
fprintf(frm, 'RMS Difference (Mean, Std, Max, Min): %3.6f, %3.6f, %3.6f, %3.6f\n',...
    rmsDiffStats.mu, rmsDiffStats.std, rmsDiffStats.max, rmsDiffStats.min);
fprintf(frm, 'RMS Percent Difference (Mean, Std, Max, Min): %3.2f, %3.2f, %3.2f, %3.2f\n\n',...
    rmsPercentDiffStats.mu, rmsPercentDiffStats.std, rmsPercentDiffStats.max, rmsPercentDiffStats.min);
fprintf(frm, 'Average FLE Model: Wishart Test Results\n-------------------------------\n');
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', passPercent3);
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', passPercent4);
fprintf(frm, 'RMS Difference (Mean, Std, Max, Min): %3.6f, %3.6f, %3.6f, %3.6f\n',...
    rmsDiffAvgFLEStats.mu, rmsDiffAvgFLEStats.std, rmsDiffAvgFLEStats.max, rmsDiffAvgFLEStats.min);
fprintf(frm, 'RMS Percent Difference (Mean, Std, Max, Min): %3.2f, %3.2f, %3.2f, %3.2f\n',...
    rmsPercentDiffAvgFLEStats.mu, rmsPercentDiffAvgFLEStats.std, rmsPercentDiffAvgFLEStats.max, rmsPercentDiffAvgFLEStats.min);


%save the parameters to be recalled later during postprocessing.
save('parm', 'parm');

% the file to screen and send an email that it is complete.
%type readme.txt;
endtime = clock;
msg = sprintf(['TRE Simulation Complete.\n    '...
    'Started:  %d-%02d-%02d %02d:%02d\n    '...
    'Ended:    %d-%02d-%02d %02d:%02d\n    '...
    'Test ID: %s\n\n    See attachment for details.'],...
    starttime(1:5), endtime(1:5), datetime);

if( bEmail )
    sendmail(email_address, 'TRE Simulation Complete', msg,...
        {'readme.txt','data.csv'} );
end

fclose('all');
cd ..\..
plotRandFiducialTargetNonHomogemousExperiment(datadir,datetime,'Non-Homogenous', 18);