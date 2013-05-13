function test = runTRESimulations(parm)

cd(parm.experimentDir);
%cd 'E:\docs\research\phd\experiments\TRE Experiments'
starttime = clock;
datetime   = sprintf('%4d%02d%02d-%02d%02d', starttime(1:5));
mkdir(datetime);
cd(datetime);
%% set up FLE components
%SigmaIso = (0.33^2)/3 * eye(3);
%SigmaAniso = diag([0.0995 0.0995 0.2985].^2);
%% perform the hypothesis test M times.
M = parm.M;
N = parm.N;

%% set up test cases.
nCases = length(parm.case);

for testcase = 1:nCases
    casename = parm.case{testcase}.name;
    Sigma = parm.case{testcase}.fleSigma;
    probe.Rigid = parm.case{testcase}.Rigid;
    
    % set up the output directory.
    mkdir(casename);
    cd(casename);

    % set up the output file.
    outFileName = sprintf('%s.csv', casename);
    outFile = fopen(outFileName,'wt');
    count = 0;

    % rotate markers into measured space.
    probe.xfrm = parm.case{testcase}.xfrm;
    
    probe.experiment.HypTestBool = zeros(M,1);
    probe.experiment.RMSPercentDiff = zeros(M,1);
    probe.experiment.PredictedRMS = zeros(M,1);
    probe.experiment.MeasuredRMS = zeros(M,1);
    while(count < M)
        if( parm.case{testcase}.varyRotation )
            screw.axis = getRandSphereSurfacePoints(1);
            screw.angle = parm.case{testcase}.varyRotationAngle * randn(1);
            probe.xfrm.R = screw2rot(screw) * probe.xfrm.R;
        end

        [testResult, probe] = simTRE(Sigma, N, probe, probe.xfrm);
        if( testResult > -1 )
            count = count + 1;
            probe.experiment.HypTestBool(count) = testResult;
            probe.experiment.PredictedRMS(count) = probe.Actual.stats.RMS;
            probe.experiment.MeasuredRMS(count) = probe.Meas.stats.RMS;
            probe.experiment.RMSPercentDiff(count) = 100*((probe.experiment.MeasuredRMS(count) - probe.experiment.PredictedRMS(count))/probe.experiment.MeasuredRMS(count));
            fprintf( outFile, '%d, %d, %3.4f, %3.4f, %3.4f\n',...
                count, probe.experiment.HypTestBool(count),...
                probe.Actual.stats.RMS, probe.Meas.stats.RMS, probe.experiment.RMSPercentDiff(count));
            fprintf('Completed %3.2f%% of test case #%d, %s\r', (100*count/M), testcase, casename);
        end
    end
        
    %% post process the result
    probe.experiment.HypTestPassPercent = 100*sum(probe.experiment.HypTestBool)/M;
    test.HypTestPassPercent(testcase) = probe.experiment.HypTestPassPercent;
    test.key(testcase) = parm.case{testcase}.key;
    fprintf( '%s: %3.1f Passed the Wishart Hypothesis Test.\n', ...
        casename, probe.experiment.HypTestPassPercent);
      
    %% set up the summary file.
    summaryFileName = sprintf( '%sSummary.txt', casename);
    summaryFile = fopen(summaryFileName, 'wt');
    fprintf(summaryFile, 'Case: %s\nDate: %s\n\n', casename, datetime);
    fprintf(summaryFile, '%3.2f%% passed the Wishart Hypothesis Test\n', probe.experiment.HypTestPassPercent);
    
    %% plot the histogram.
    figure(1);
    hist(probe.experiment.RMSPercentDiff);
    title('Histogram of the Percentage Differences for RMS', 'fontsize',14);
    set(gca,'fontsize',14);
    figurefilename = sprintf('Histogram_%s', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);
    
    %% plot the predicted v. measured plot.
    figure(2);
    plot(probe.experiment.PredictedRMS, probe.experiment.MeasuredRMS, '.');
    hold on;
    probe.experiment.meanRMSPercentDiff = mean(probe.experiment.RMSPercentDiff);
    probe.experiment.stdRMSPercentDiff = std(probe.experiment.RMSPercentDiff);
    probe.experiment.maxRMSPercentDiff = max(probe.experiment.RMSPercentDiff);
    probe.experiment.minRMSPercentDiff = min(probe.experiment.RMSPercentDiff);
    fprintf(summaryFile, 'RMS Percent Difference:\n Mean:\t%3.2f%%\n Std Dev:%3.2f%%\t\n Max:\t%3.2f%%\n Min:\t%3.2f%%\n',...
        probe.experiment.meanRMSPercentDiff, probe.experiment.stdRMSPercentDiff,...
        probe.experiment.maxRMSPercentDiff, probe.experiment.minRMSPercentDiff);
    maxValue = max([probe.experiment.PredictedRMS; probe.experiment.MeasuredRMS]);
    minValue = min([probe.experiment.PredictedRMS; probe.experiment.MeasuredRMS]);
    plot([minValue maxValue], [minValue, maxValue], 'k');
    xlabel('Predicted TRE RMS (mm)', 'fontsize',14);
    ylabel('Simulated TRE RMS (mm)', 'fontsize',14);
    set(gca,'fontsize',14);
    title('Comparison of Predicted and Simulated TRE RMS', 'fontsize',14);
    hold off;
    axis([minValue, maxValue, minValue, maxValue]);
    figurefilename = sprintf('PredvMeas_%s', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);

    fclose(outFile);
    fclose(summaryFile);
    %% save the data.
    save(casename, 'probe');
    cd ..
    figure(3);
    plot(test.key, test.HypTestPassPercent);
end

save(datetime, 'parm', 'test');
fclose('all');
cd ..

figure(3);
plot(test.key, test.HypTestPassPercent);