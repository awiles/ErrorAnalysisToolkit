function [probe,rmsDiff] = runRandFiducialTargetExperiment(nSamples, nBodies, nTrials, nOrientations, nPositions,...
    nMarkers, Sigma, mrkRange, targetRange, varargin)
% [probe,rmsDiff] = runRandFiducialTargetExperiment(nSamples, nTrials,...
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
%       Sigma    - FLE given as a covariance matrix.
%       mrkRange - distance over which the markers will appear.  The origin
%                       is assumed to bisect this range.
%       targetRange - distance over which the target will appear.  The
%                       origin is assumed to bisect this range.

% defaults for optional arguments.
smtp_server = '';
email_address = '';
bEmail = 0;

if( nargin > 8 )
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
    bEMail = 1;
end

nCount=0;
nTotalCount = nBodies*nTrials*nOrientations*nPositions;
passCount1 = 0;
passCount2 = 0;
rmsDiff = zeros(nTotalCount,3);

% set up the data directory to store the results.
starttime = clock;
datetime   = sprintf('%4d%02d%02d-%02d%02d', starttime(1:5));
mkdir(datetime);
cd(datetime);

% write the experiment parameters to the readme.txt file.
frm = fopen('readme.txt', 'wt');
fprintf(frm, 'Date processed: %s\n', datetime); parm.name = datetime;
fprintf(frm, 'nSamples = %d\n', nSamples); parm.nSamples = nSamples;
fprintf(frm, 'nBodies = %d\n', nBodies); parm.nBodies = nBodies;
fprintf(frm, 'nTrials = %d\n', nTrials); parm.nTrials = nTrials;
fprintf(frm, 'nOrientations = %d\n', nOrientations); parm.nOrientations = nOrientations;
fprintf(frm, 'nPositions = %d\n', nPositions); parm.nPositions = nPositions;
fprintf(frm, 'nMarkers = %d\n', nMarkers); parm.nMarkers = nMarkers;
fprintf(frm, 'mrkRange = %d\n', mrkRange); parm.mrkRange = mrkRange;
fprintf(frm, 'targetRange = %d\n\n', targetRange); parm.targetRange = targetRange;
fprintf(frm, 'Sigma = \n');
fprintf(frm, '        %f %f %f \n', Sigma); parm.Sigma = Sigma;

fprintf(frm, 'Data Summary\n------------\n\n');

% set up the data file.
fdata = fopen('data.csv', 'wt');

% perform the simulation.
for i=1:nBodies
    % build the marker set.
    probe = generateRandFiducialTargets(nMarkers, mrkRange, targetRange);
    probe.bHomogenous=1;
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
                %assuming origin for now.
                %probeXfrm.R=eye(3);
                probeXfrm.pos = [0 0 1000];
                %****************Transform into the Random Xfrm*******%
                %move probe coordinate frame to tip.
                %temp.Rigid.mrk = probe.Rigid.mrk - repmat(probe.Rigid.tip, size(probe.Rigid.mrk,1), 1);
                %temp.Rigid.tip = [0 0 0];
                %transform probe rigid body into test space.
                probe.Actual.mrk = (probeXfrm.R * probe.Rigid.mrk')'...
                    + repmat(probeXfrm.pos, size(probe.Rigid.mrk,1), 1);
                probe.Actual.tip = (probeXfrm.R * probe.Rigid.tip')' + probeXfrm.pos;

                %**************Compute the theoretical statistics.*******%
                probe.Actual.stats.mu = zeros(1,3);
                [probe.Actual.stats.RMS, probe.Actual.stats.Sigma, probe.Actual.stats.SigmaPA] =...
                    calcTRE(Sigma, [probe.Actual.mrk; probe.Actual.tip]);
                probe.Fitz.stats.RMS = calcTRE(sqrt(trace(Sigma)), [probe.Actual.mrk; probe.Actual.tip]);

                %*******************Monte Carlo trial.**********%
                [probe.Meas.error, probe.Meas.tip] = simTRE(Sigma, nSamples, probe.Rigid, probe.Actual );
                probe = computeStats(probe);
                %compare the simulation to the theoretical stats.
                testResult = compareTREStats(probe);
                passCount1 = passCount1 + testResult.probe.covariance;
                passCount2 = passCount2 + testResult.probe.meanandcov;
                nCount = nCount + 1;
                rmsDiff(nCount,1) = testResult.probe.RMS.Meas;
                rmsDiff(nCount,2) = testResult.probe.RMS.Theory;
                rmsDiff(nCount,3) = testResult.probe.RMS.Diff;
                rmsDiff(nCount,4) = testResult.probe.RMS.PercentDiff;                
                fprintf('Completed %3.1f%% ... nMarkers = %d \n', (100*nCount/nTotalCount), nMarkers);
                
                % save the data.
                filename = sprintf('data%06d', nCount);
                save(filename, 'probe', 'testResult');
                % write out the summary data.
                fprintf(fdata,'%d, %d, %d, %3.6f, %3.6f, %3.6f, %3.2f\n', nCount, ...
                    testResult.probe.covariance, testResult.probe.meanandcov,...
                    testResult.probe.RMS.Meas, testResult.probe.RMS.Theory,...
                    testResult.probe.RMS.Diff, testResult.probe.RMS.PercentDiff);
            end
        end
    end
end

%compute the percentage of successful Wishart tests.
passPercent1(i) = 100*passCount1/nCount;
passPercent2(i) = 100*passCount2/nCount;
% compute the stats of the RMS difference.
rmsDiffStats.mu = mean(rmsDiff(:,3));
rmsDiffStats.std = std(rmsDiff(:,3));
rmsDiffStats.max = max(rmsDiff(:,3));
rmsDiffStats.min = min(rmsDiff(:,3));
% compute the stats of the RMS percent difference.
rmsPercentDiffStats.mu = mean(rmsDiff(:,4));
rmsPercentDiffStats.std = std(rmsDiff(:,4));
rmsPercentDiffStats.max = max(rmsDiff(:,4));
rmsPercentDiffStats.min = min(rmsDiff(:,4));

% write summary to readme file.
fprintf(frm, 'Wishart Test Results, Trial #%d\n-------------------------------\n', i);
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', passPercent1(i));
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', passPercent2(i));
fprintf(frm, 'RMS Difference (Mean, Std, Max, Min): %3.6f, %3.6f, %3.6f, %3.6f\n',...
    rmsDiffStats.mu, rmsDiffStats.std, rmsDiffStats.max, rmsDiffStats.min);
fprintf(frm, 'RMS Percent Difference (Mean, Std, Max, Min): %3.2f, %3.2f, %3.2f, %3.2f\n',...
    rmsPercentDiffStats.mu, rmsPercentDiffStats.std, rmsPercentDiffStats.max, rmsPercentDiffStats.min);
%save the parameters to be recalled later during postprocessing.
save('parm', 'parm');

% the file to screen and send an email that it is complete.
type readme.txt;
endtime = clock;
msg = sprintf(['TRE Simulation Complete.\n    '...
    'Started:  %d-%02d-%02d %02d:%02d\n    '...
    'Ended:    %d-%02d-%02d %02d:%02d\n    '...
    'Test ID: %s\n\n    See attachment for details.'],...
    starttime(1:5), endtime(1:5), datetime);

if( bEmail)
    sendmail(email_address, 'TRE Simulation Complete', msg, {'readme.txt','data.csv'});
end

fclose('all');
cd ..
plotRandFiducialTargetExperiment(datetime,'Homogenous', 18);