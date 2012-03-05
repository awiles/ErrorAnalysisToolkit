function [probe,rmsDiff] = runRandProbeRefXfrmExperiment(nSamples, nOrientations, nPositions,...
    design, A, B, rho, r, d, Sigma)
% [probe,rmsDiff] = runRandProbeXfrmRefExperiment(nSamples, nOrientations, nPositions,...
%                       design, A, B, rho, r, d, Sigma)
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
passCount1 = 0; passCount2 = 0; passCount3 = 0;
passCount4 = 0; passCount5 = 0; passCount6 = 0;
rmsDiff = zeros(nTotalCount,12);

% set up the data directory to store the results.
starttime = clock;
datetime   = sprintf('%4d%02d%02d-%02d%02d', starttime(1:5));
mkdir(datetime);
cd(datetime);

% write the experiment parameters to the readme.txt file.
frm = fopen('readme.txt', 'wt');
fprintf(frm, 'Random Xfrm Experiment using a Reference\n');
fprintf(frm, 'Date processed: %s\n', datetime); parm.name = datetime;
fprintf(frm, 'nSamples = %d\n', nSamples); parm.nSamples = nSamples;
fprintf(frm, 'nOrientations = %d\n', nOrientations); parm.nOrientations = nOrientations;
fprintf(frm, 'nPositions = %d\n', nPositions); parm.nPositions = nPositions;
fprintf(frm, 'West Design = %s\n', design); parm.design = design;
fprintf(frm, 'A = %3.2f\n', A); parm.A = A;
fprintf(frm, 'B = %3.2f\n', B); parm.B = B;
fprintf(frm, 'rho = %3.2f\n', rho); parm.rho = rho;
fprintf(frm, 'r = %3.2f\n', r); parm.r = r;
fprintf(frm, 'd = %3.2f\n', d); parm.d = d;
fprintf(frm, 'Sigma = \n');
fprintf(frm, '        %f %f %f \n', Sigma); parm.Sigma = Sigma;

fprintf(frm, 'Data Summary\n------------\n\n');

% set up the data file.
fdata = fopen('data.csv', 'wt');

% perform the simulation.
% build the marker set.
[probe.Rigid.mrk, probe.Rigid.normals, probe.Rigid.tip] = getWestToolDesign(design, A, B, rho);
[ref.Rigid.mrk, ref.Rigid.normals, ref.Rigid.tip] = getWestToolDesign('e', r, r, 0);

for k=1:nOrientations
    for m=1:nPositions
        %******************Get the Random Transform ***********%
        probe = generateRandProbeXfrm(probe, 60);
        ref = generateRandProbeXfrm(ref, 60);

        %*********Move Reference to a Position d mm from tool tip ****%
        goodRef = 0;
        while(~goodRef)
            randDirVec = rand(1,3);
            randDirVec = (1/(sqrt(randDirVec*randDirVec'))) * randDirVec;
            if( (acosd(randDirVec(1,3)) >= 60) & (acosd(randDirVec(1,3)) <= 120) )
                goodRef = 1;
            end
        end
        ref.xfrm.pos = probe.Actual.tip + d*randDirVec;

        % transform rigid body into test space.
        ref.Actual.mrk = (ref.xfrm.R * ref.Rigid.mrk')'...
            + repmat(ref.xfrm.pos, size(ref.Rigid.mrk,1), 1);
        ref.Actual.normals = (ref.xfrm.R * ref.Rigid.normals')';
        ref.Actual.tip = probe.Actual.tip;
        ref.Rigid.tip = (ref.xfrm.R' * (ref.Actual.tip - ref.xfrm.pos)')';

        %DEBUG:
        %dMeas = sqrt(ref.Rigid.tip*ref.Rigid.tip');
        %fprintf('DEBUG: distance to target in ref space: %3.2f\n', dMeas);

        %**************Compute the theoretical statistics.*******%
        probe.Actual.stats.mu = zeros(1,3);
        [probe.Actual.stats.RMS, probe.Actual.stats.Sigma, probe.Actual.stats.SigmaPA] =...
            calcTRE(Sigma, [probe.Actual.mrk; probe.Actual.tip], 0);
        probe.Fitz.stats.RMS = calcTREold(sqrt(trace(Sigma)), [probe.Actual.mrk; probe.Actual.tip]);

        ref.Actual.stats.mu = zeros(1,3);
        [ref.Actual.stats.RMS, ref.Actual.stats.Sigma, ref.Actual.stats.SigmaPA] =...
            calcTRE(Sigma, [ref.Actual.mrk; ref.Actual.tip], 0);
        ref.Fitz.stats.RMS = calcTREold(sqrt(trace(Sigma)), [ref.Actual.mrk; ref.Actual.tip]);

        combined.Actual.stats.mu = probe.Actual.stats.mu + ref.Actual.stats.mu;
        combined.Actual.stats.RMS = sqrt(probe.Actual.stats.RMS^2 + ref.Actual.stats.RMS^2);
        combined.Actual.stats.Sigma = ref.xfrm.R'*(ref.Actual.stats.Sigma + probe.Actual.stats.Sigma)*ref.xfrm.R;

        %*******************Monte Carlo trial.**********%
        [probe.Meas.error, probe.Meas.tip, ref.Meas.error, ref.Meas.tip, ...
            combined.Meas.error, combined.Meas.tip] =...
            simTRE(Sigma, nSamples, probe.Rigid, probe.Actual, ref.Rigid, ref.Actual );
        probe = computeStats(probe);
        ref = computeStats(ref);
        combined = computeStats(combined);
        
        %compare the simulation to the theoretical stats.
        testResult = compareTREStats(probe, ref, combined);
        passCount1 = passCount1 + testResult.probe.covariance;
        passCount2 = passCount2 + testResult.probe.meanandcov;
        passCount3 = passCount3 + testResult.ref.covariance;
        passCount4 = passCount4 + testResult.ref.meanandcov;
        passCount5 = passCount5 + testResult.combined.covariance;
        passCount6 = passCount6 + testResult.combined.meanandcov;
        nCount = nCount + 1;
        
        rmsDiff(nCount,1) = testResult.probe.RMS.Meas;
        rmsDiff(nCount,2) = testResult.probe.RMS.Theory;
        rmsDiff(nCount,3) = testResult.probe.RMS.Diff;
        rmsDiff(nCount,4) = testResult.probe.RMS.PercentDiff;
        rmsDiff(nCount,5) = testResult.ref.RMS.Meas;
        rmsDiff(nCount,6) = testResult.ref.RMS.Theory;
        rmsDiff(nCount,7) = testResult.ref.RMS.Diff;
        rmsDiff(nCount,8) = testResult.ref.RMS.PercentDiff;
        rmsDiff(nCount,9) = testResult.combined.RMS.Meas;
        rmsDiff(nCount,10) = testResult.combined.RMS.Theory;
        rmsDiff(nCount,11) = testResult.combined.RMS.Diff;
        rmsDiff(nCount,12) = testResult.combined.RMS.PercentDiff;
        fprintf('Completed %3.1f%% ... design = %s, rho = %3.2f, r = %3.2f, d = %3.2f \n',...
            (100*nCount/nTotalCount), design, rho, r, d);

        % save the data.
        filename = sprintf('data%06d', nCount);
        save(filename, 'probe', 'ref', 'combined','testResult');
        % write out the summary data.
        fprintf(fdata,'%d, %d, %d, %3.6f, %3.6f, %3.6f, %3.2f, %d, %d, %3.6f, %3.6f, %3.6f, %3.2f, %d, %d, %3.6f, %3.6f, %3.6f, %3.2f \n', nCount, ...
            testResult.probe.covariance, testResult.probe.meanandcov,...
            testResult.probe.RMS.Meas, testResult.probe.RMS.Theory,...
            testResult.probe.RMS.Diff, testResult.probe.RMS.PercentDiff,...
            testResult.ref.covariance, testResult.ref.meanandcov,...
            testResult.ref.RMS.Meas, testResult.ref.RMS.Theory,...
            testResult.ref.RMS.Diff, testResult.ref.RMS.PercentDiff,...
            testResult.combined.covariance, testResult.combined.meanandcov,...
            testResult.combined.RMS.Meas, testResult.combined.RMS.Theory,...
            testResult.combined.RMS.Diff, testResult.combined.RMS.PercentDiff);
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
fprintf(frm, '\nProbe Results\n');
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', 100*passCount1/nCount);
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', 100*passCount2/nCount);
fprintf(frm, 'RMS Difference (Mean, Std, Max, Min): %3.6f, %3.6f, %3.6f, %3.6f\n',...
    mean(rmsDiff(:,3)), std(rmsDiff(:,3)), max(rmsDiff(:,3)), min(rmsDiff(:,3)));
fprintf(frm, 'RMS Percent Difference (Mean, Std, Max, Min): %3.2f, %3.2f, %3.2f, %3.2f\n',...
    mean(rmsDiff(:,4)), std(rmsDiff(:,4)), max(rmsDiff(:,4)), min(rmsDiff(:,4)));
fprintf(frm, '\nReference Results\n');
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', 100*passCount3/nCount);
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', 100*passCount4/nCount);
fprintf(frm, 'RMS Difference (Mean, Std, Max, Min): %3.6f, %3.6f, %3.6f, %3.6f\n',...
    mean(rmsDiff(:,7)), std(rmsDiff(:,7)), max(rmsDiff(:,7)), min(rmsDiff(:,7)));
fprintf(frm, 'RMS Percent Difference (Mean, Std, Max, Min): %3.2f, %3.2f, %3.2f, %3.2f\n',...
    mean(rmsDiff(:,8)), std(rmsDiff(:,8)), max(rmsDiff(:,8)), min(rmsDiff(:,8)));
fprintf(frm, '\nCombined Results\n');
fprintf(frm, 'Covariance Matrix Only: %3.2f\n', 100*passCount5/nCount);
fprintf(frm, 'Mean and Covariance Matrix: %3.2f\n', 100*passCount6/nCount);
fprintf(frm, 'RMS Difference (Mean, Std, Max, Min): %3.6f, %3.6f, %3.6f, %3.6f\n',...
    mean(rmsDiff(:,11)), std(rmsDiff(:,11)), max(rmsDiff(:,11)), min(rmsDiff(:,11)));
fprintf(frm, 'RMS Percent Difference (Mean, Std, Max, Min): %3.2f, %3.2f, %3.2f, %3.2f\n',...
    mean(rmsDiff(:,12)), std(rmsDiff(:,12)), max(rmsDiff(:,12)), min(rmsDiff(:,12)));


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