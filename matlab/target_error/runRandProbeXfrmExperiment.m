function [probe,rmsDiff] = runRandProbeXfrmExperiment(nSamples, nOrientations, nPositions,...
    design, A, B, rho, Sigma)
% [probe,rmsDiff] = runRandProbeXfrmExperiment(nSamples, nOrientations, nPositions,...
%                       design, A, B, rho, Sigma)
%
%       nSamples - number of random samples in each Monte Carlo simulation,
%                   i.e. N = 10,000 or 10^6 or something similar.
%       nOrientations - number of random orientations to test for each
%                   trial.
%       nPositions - number of random positions to test (N/A if Sigma is
%                       constant throughout the volume, therefore set to
%                       one.  Kept here as a place holder for future work.)
%       Sigma    - FLE given as a covariance matrix.

email_setup;

nCount=0;
nTotalCount = nOrientations*nPositions;
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
fprintf(frm, 'Random Xfrm Experiment\n');
fprintf(frm, 'Date processed: %s\n', datetime); parm.name = datetime;
fprintf(frm, 'nSamples = %d\n', nSamples); parm.nSamples = nSamples;
fprintf(frm, 'nOrientations = %d\n', nOrientations); parm.nOrientations = nOrientations;
fprintf(frm, 'nPositions = %d\n', nPositions); parm.nPositions = nPositions;
fprintf(frm, 'West Design = %s\n', design); parm.design = design;
fprintf(frm, 'A = %3.2f\n', A); parm.A = A;
fprintf(frm, 'B = %3.2f\n', B); parm.B = B;
fprintf(frm, 'rho = %3.2f\n', rho); parm.rho = rho;
fprintf(frm, 'Sigma = \n');
fprintf(frm, '        %f %f %f \n', Sigma); parm.Sigma = Sigma;

fprintf(frm, 'Data Summary\n------------\n\n');

% set up the data file.
fdata = fopen('data.csv', 'wt');

% perform the simulation.
% build the marker set.
[probe.Rigid.mrk, probe.Rigid.normals, probe.Rigid.tip] = getWestToolDesign(design, A, B, rho);

for k=1:nOrientations    
    for m=1:nPositions
        %******************Get the Random Transform ***********%
        probe = generateRandProbeXfrm(probe, 60);
        
        %**************Compute the theoretical statistics.*******%
        probe.Actual.stats.mu = zeros(1,3);
        [probe.Actual.stats.RMS, probe.Actual.stats.Sigma, probe.Actual.stats.SigmaPA] =...
            calcTRE(Sigma, [probe.Actual.mrk; probe.Actual.tip], 0);
        probe.Fitz.stats.RMS = calcTREold(sqrt(trace(Sigma)), [probe.Actual.mrk; probe.Actual.tip]);

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
        fprintf('Completed %3.1f%% ... design = %s, rho = %3.2f \n', (100*nCount/nTotalCount), design, rho);

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

%compute the percentage of successful Wishart tests.
passPercent1 = 100*passCount1/nCount;
passPercent2 = 100*passCount2/nCount;
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
fprintf(frm, 'Wishart Test Results\n-------------------------------\n');
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', passPercent1);
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', passPercent2);
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

sendmail('awiles@imaging.robarts.ca', 'TRE Simulation Complete', msg, {'readme.txt','data.csv'});
fclose('all');
cd ..