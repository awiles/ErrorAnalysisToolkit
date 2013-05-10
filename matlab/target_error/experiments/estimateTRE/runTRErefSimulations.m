clear all;
close all;

cd 'E:\docs\research\phd\experiments\TRE Experiments'
starttime = clock;
datetime   = sprintf('%4d%02d%02d-%02d%02d', starttime(1:5));
mkdir(datetime);
cd(datetime);
%% set up FLE components
SigmaIso = (0.33^2)/3 * eye(3);
SigmaAniso = diag([0.0995 0.0995 0.2985].^2);
%% perform the hypothesis test M times.
M = 1000;
N = 10000;

%% set up probe marker coordinates.  Rigid body parameters.
% define West tool design parameters.


for testcase = 4
    switch(testcase)
        case 1
            casename = 'Des_e_Iso_85mmTip_d100_r32';
            probeDesign = 'e';
            Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 100;
        case 2
            casename = 'Des_e_Iso_85mmTip_d200_r32';
            probeDesign = 'e';
            Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 200;
        case 3
            casename = 'Des_e_Iso_85mmTip_d300_r32';
            probeDesign = 'e';
            Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 300;
        case 4
            casename = 'Des_e_Iso_85mmTip_d400_r32';
            probeDesign = 'e';
            Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 400;
        case 5
            casename = 'Des_e_Aniso_85mmTip_d100_r32';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 100;
        case 6
            casename = 'Des_e_Aniso_85mmTip_d200_r32';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 200;
        case 7
            casename = 'Des_e_Aniso_85mmTip_d300_r32';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 300;
        case 8
            casename = 'Des_e_Aniso_85mmTip_d400_r32';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
            d = 400;
        case 9
            casename = 'Des_e_Aniso_85mmTip_d100_r64';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 64;
            d = 100;
        case 10
            casename = 'Des_e_Aniso_85mmTip_d200_r64';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 64;
            d = 200;
        case 11
            casename = 'Des_e_Aniso_85mmTip_d300_r64';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 64;
            d = 300;
        case 12
            casename = 'Des_e_Aniso_85mmTip_d400_r64';
            probeDesign = 'e';
            Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 64;
            d = 400;
        otherwise
            error('Ummm... that case does not exist');
    end

    mkdir(casename);
    cd(casename);

    switch('d')
        case {'d'}
            % West design in Fig. 2(d).
            probe.Rigid.mrk = [ 0 B/2 0;
                -A/2 0 0;
                0 -B/2 0;
                A/2 0 0 ];
            probe.Rigid.normals = repmat([ 0 0 1], 4, 1);
        case{'e'}
            % West design in Fig. 2(e).
            probe.Rigid.mrk = [ -A/2 B/2 0;
                -A/2 -B/2 0;
                A/2 -B/2 0;
                A/2 B/2 0 ];
            probe.Rigid.normals = repmat([ 0 0 1], 4, 1);
        otherwise
            error('Invalid tool design case.');
    end
    % tip is located at
    probe.Rigid.tip = [ -rho 0 0 ];

    %% set up reference marker coordinates.  Rigid body parameters.
    % To be used in later experiments.

    % West design in Fig. 2(e).
    ref.Rigid.mrk = [ -r/2 r/2 0;
        -r/2 -r/2 0;
        r/2 -r/2 0;
        r/2 r/2 0 ];
    ref.Rigid.normals = repmat([ 0 0 1], 4, 1);
    % tip is located at
    ref.Rigid.tip = [ 0 0 0 ];

    outFileName = sprintf('%s.csv', casename);
    outFile = fopen(outFileName,'wt');
    count = 0;

    % rotate markers into measured space.
    probePos = [0 100 -1000];
    rotX = -10; rotY = 0; rotZ = 10;
    % R = getRotMatrixd([rotX, rotY, rotZ]);
    refPos = [0, -200, -1000];
    refR = getRotMatrixd([rotX, rotY, rotZ]);
    testBool = zeros(M,6);
    
    probeRMSPercentDiff = zeros(M,1);
    probePredictedRMS = zeros(M,1);
    probeMeasuredRMS = zeros(M,1);
    probeXSkewness = zeros(M,1);
    probeYSkewness = zeros(M,1);
    probeZSkewness = zeros(M,1);
    probeXKurtosis = zeros(M,1);
    probeYKurtosis = zeros(M,1);
    probeZKurtosis = zeros(M,1);
    
    refRMSPercentDiff = zeros(M,1);
    refPredictedRMS = zeros(M,1);
    refMeasuredRMS = zeros(M,1);
    refXSkewness = zeros(M,1);
    refYSkewness = zeros(M,1);
    refZSkewness = zeros(M,1);
    refXKurtosis = zeros(M,1);
    refYKurtosis = zeros(M,1);
    refZKurtosis = zeros(M,1);
    
    combinedRMSPercentDiff = zeros(M,1);
    combinedPredictedRMS = zeros(M,1);
    combinedMeasuredRMS = zeros(M,1);
    combinedXSkewness = zeros(M,1);
    combinedYSkewness = zeros(M,1);
    combinedZSkewness = zeros(M,1);
    combinedXKurtosis = zeros(M,1);
    combinedYKurtosis = zeros(M,1);
    combinedZKurtosis = zeros(M,1);
    
    while(count < M)
        % get random probe position and orientation.
        screw.axis = getRandSphereSurfacePoints(1);
        screw.angle = 45 * randn(1);
        probeXfrm.R = screw2rot(screw);
        probeXfrm.pos = probePos;
        % get random reference position and orientation.
        goodReference = 0;
        while( ~goodReference )
            screw.axis = getRandSphereSurfacePoints(1);
            screw.angle = 15 * randn(1);
            refXfrm.R = screw2rot(screw);
            randDirVec = rand(1,3);
            randDirVec = (1/(sqrt(randDirVec*randDirVec'))) * randDirVec;
            if( (acosd(randDirVec(1,3)) >= 60) & (acosd(randDirVec(1,3)) <= 120) )
                goodReference = 1;
                break;
            end
        end
        refXfrm.pos = probePos + d*randDirVec;
        [testResult, probe, ref, combined] = simTREref(Sigma, N, probe, probeXfrm, ref, refXfrm);
        if( testResult.probe.covariance > -1 )
            count = count + 1;
            testBool(count,:) = [testResult.probe.covariance, testResult.ref.covariance, testResult.combined.covariance, ...
                testResult.probe.meanandcov, testResult.ref.meanandcov, testResult.combined.meanandcov];
            
            probePredictedRMS(count) = probe.Actual.stats.RMS;
            probeMeasuredRMS(count) = probe.Meas.stats.RMS;
            probeRMSPercentDiff(count)...
                = 100*((probeMeasuredRMS(count) - probePredictedRMS(count))/probeMeasuredRMS(count));
            probeXSkewness(count) = probe.Meas.stats.skewness.x;
            probeYSkewness(count) = probe.Meas.stats.skewness.y;
            probeZSkewness(count) = probe.Meas.stats.skewness.z;
            probeXKurtosis(count) = probe.Meas.stats.kurtosis.x;
            probeYKurtosis(count) = probe.Meas.stats.kurtosis.y;
            probeZKurtosis(count) = probe.Meas.stats.kurtosis.z;
            
            refPredictedRMS(count) = ref.Actual.stats.RMS;
            refMeasuredRMS(count) = ref.Meas.stats.RMS;
            refRMSPercentDiff(count)...
                = 100*((refMeasuredRMS(count) - refPredictedRMS(count))/refMeasuredRMS(count));
            refXSkewness(count) = ref.Meas.stats.skewness.x;
            refYSkewness(count) = ref.Meas.stats.skewness.y;
            refZSkewness(count) = ref.Meas.stats.skewness.z;
            refXKurtosis(count) = ref.Meas.stats.kurtosis.x;
            refYKurtosis(count) = ref.Meas.stats.kurtosis.y;
            refZKurtosis(count) = ref.Meas.stats.kurtosis.z;
            
            combinedPredictedRMS(count) = combined.Actual.stats.RMS;
            combinedMeasuredRMS(count) = combined.Meas.stats.RMS;
            combinedRMSPercentDiff(count)...
                = 100*((combinedMeasuredRMS(count) - combinedPredictedRMS(count))/combinedMeasuredRMS(count));
            combinedXSkewness(count) = combined.Meas.stats.skewness.x;
            combinedYSkewness(count) = combined.Meas.stats.skewness.y;
            combinedZSkewness(count) = combined.Meas.stats.skewness.z;
            combinedXKurtosis(count) = combined.Meas.stats.kurtosis.x;
            combinedYKurtosis(count) = combined.Meas.stats.kurtosis.y;
            combinedZKurtosis(count) = combined.Meas.stats.kurtosis.z;
            
            fprintf( outFile, '%d, %d, %d, %d, %d, %d, %d, %3.4f, %3.4f, %3.4f, %3.4f, %3.4f, %3.4f, %3.4f, %3.4f, %3.4f\n',...
                count, testResult.probe.covariance, testResult.ref.covariance, testResult.combined.covariance,...
                testResult.probe.meanandcov, testResult.ref.meanandcov, testResult.combined.meanandcov, ...
                probe.Actual.stats.RMS, probe.Meas.stats.RMS, probeRMSPercentDiff(count),...
                ref.Actual.stats.RMS, ref.Meas.stats.RMS, refRMSPercentDiff(count),...
                combined.Actual.stats.RMS, combined.Meas.stats.RMS, combinedRMSPercentDiff(count));
            fprintf('Completed %3.2f%% of test case #%d, %s\r', (100*count/M), testcase, casename);
        end
    end

    passPercent = 100/M*sum(testBool)
    summaryFile = fopen('summary.txt', 'wt');
    fprintf(summaryFile, 'Probe    : %3.2f%% passed the Wishart Hypothesis Test (Covariance Only)\n', passPercent(1));
    fprintf(summaryFile, 'Reference: %3.2f%% passed the Wishart Hypothesis Test (Covariance Only)\n', passPercent(2));
    fprintf(summaryFile, 'Combined : %3.2f%% passed the Wishart Hypothesis Test (Covariance Only)\n', passPercent(3));
    fprintf(summaryFile, 'Probe    : %3.2f%% passed the Wishart Hypothesis Test (Mean and Covariance)\n', passPercent(4));
    fprintf(summaryFile, 'Reference: %3.2f%% passed the Wishart Hypothesis Test (Mean and Covariance)\n', passPercent(5));
    fprintf(summaryFile, 'Combined : %3.2f%% passed the Wishart Hypothesis Test (Mean and Covariance)\n', passPercent(6));

    figure(1);
    hist(probeRMSPercentDiff);
    title('Histogram of the Percentage Differences for RMS of Probe', 'fontsize',14);
    set(gca,'fontsize',14);
    figurefilename = sprintf('Histogram_%s_Probe', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);

    figure(2);
    hist(refRMSPercentDiff);
    title('Histogram of the Percentage Differences for RMS of Reference', 'fontsize',14);
    set(gca,'fontsize',14);
    figurefilename = sprintf('Histogram_%s_Reference', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);
    
    figure(3);
    hist(combinedRMSPercentDiff);
    title('Histogram of the Percentage Differences for RMS Combined', 'fontsize',14);
    set(gca,'fontsize',14);
    figurefilename = sprintf('Histogram_%s', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);
    
    meanProbeRMSPercentDiff = mean(probeRMSPercentDiff);
    stdProbeRMSPercentDiff = std(probeRMSPercentDiff);
    maxProbeRMSPercentDiff = max(probeRMSPercentDiff);
    minProbeRMSPercentDiff = min(probeRMSPercentDiff);
    meanProbeXSkewness = mean(probeXSkewness);
    stdProbeXSkewness = std(probeXSkewness);
    maxProbeXSkewness = max(probeXSkewness);
    minProbeXSkewness = min(probeXSkewness);
    meanProbeYSkewness = mean(probeYSkewness);
    stdProbeYSkewness = std(probeYSkewness);
    maxProbeYSkewness = max(probeYSkewness);
    minProbeYSkewness = min(probeYSkewness);
    meanProbeZSkewness = mean(probeZSkewness);
    stdProbeZSkewness = std(probeZSkewness);
    maxProbeZSkewness = max(probeZSkewness);
    minProbeZSkewness = min(probeZSkewness);
    meanProbeXKurtosis = mean(probeXKurtosis);
    stdProbeXKurtosis = std(probeXKurtosis);
    maxProbeXKurtosis = max(probeXKurtosis);
    minProbeXKurtosis = min(probeXKurtosis);
    meanProbeYKurtosis = mean(probeYKurtosis);
    stdProbeYKurtosis = std(probeYKurtosis);
    maxProbeYKurtosis = max(probeYKurtosis);
    minProbeYKurtosis = min(probeYKurtosis);
    meanProbeZKurtosis = mean(probeZKurtosis);
    stdProbeZKurtosis = std(probeZKurtosis);
    maxProbeZKurtosis = max(probeZKurtosis);
    minProbeZKurtosis = min(probeZKurtosis);
    
    meanRefRMSPercentDiff = mean(refRMSPercentDiff);
    stdRefRMSPercentDiff = std(refRMSPercentDiff);
    maxRefRMSPercentDiff = max(refRMSPercentDiff);
    minRefRMSPercentDiff = min(refRMSPercentDiff);
    meanRefXSkewness = mean(refXSkewness);
    stdRefXSkewness = std(refXSkewness);
    maxRefXSkewness = max(refXSkewness);
    minRefXSkewness = min(refXSkewness);
    meanRefYSkewness = mean(refYSkewness);
    stdRefYSkewness = std(refYSkewness);
    maxRefYSkewness = max(refYSkewness);
    minRefYSkewness = min(refYSkewness);
    meanRefZSkewness = mean(refZSkewness);
    stdRefZSkewness = std(refZSkewness);
    maxRefZSkewness = max(refZSkewness);
    minRefZSkewness = min(refZSkewness);
    meanRefXKurtosis = mean(refXKurtosis);
    stdRefXKurtosis = std(refXKurtosis);
    maxRefXKurtosis = max(refXKurtosis);
    minRefXKurtosis = min(refXKurtosis);
    meanRefYKurtosis = mean(refYKurtosis);
    stdRefYKurtosis = std(refYKurtosis);
    maxRefYKurtosis = max(refYKurtosis);
    minRefYKurtosis = min(refYKurtosis);
    meanRefZKurtosis = mean(refZKurtosis);
    stdRefZKurtosis = std(refZKurtosis);
    maxRefZKurtosis = max(refZKurtosis);
    minRefZKurtosis = min(refZKurtosis);
    
    meanCombinedRMSPercentDiff = mean(combinedRMSPercentDiff);
    stdCombinedRMSPercentDiff = std(combinedRMSPercentDiff);
    maxCombinedRMSPercentDiff = max(combinedRMSPercentDiff);
    minCombinedRMSPercentDiff = min(combinedRMSPercentDiff);
    meanCombinedXSkewness = mean(combinedXSkewness);
    stdCombinedXSkewness = std(combinedXSkewness);
    maxCombinedXSkewness = max(combinedXSkewness);
    minCombinedXSkewness = min(combinedXSkewness);
    meanCombinedYSkewness = mean(combinedYSkewness);
    stdCombinedYSkewness = std(combinedYSkewness);
    maxCombinedYSkewness = max(combinedYSkewness);
    minCombinedYSkewness = min(combinedYSkewness);
    meanCombinedZSkewness = mean(combinedZSkewness);
    stdCombinedZSkewness = std(combinedZSkewness);
    maxCombinedZSkewness = max(combinedZSkewness);
    minCombinedZSkewness = min(combinedZSkewness);
    meanCombinedXKurtosis = mean(combinedXKurtosis);
    stdCombinedXKurtosis = std(combinedXKurtosis);
    maxCombinedXKurtosis = max(combinedXKurtosis);
    minCombinedXKurtosis = min(combinedXKurtosis);
    meanCombinedYKurtosis = mean(combinedYKurtosis);
    stdCombinedYKurtosis = std(combinedYKurtosis);
    maxCombinedYKurtosis = max(combinedYKurtosis);
    minCombinedYKurtosis = min(combinedYKurtosis);
    meanCombinedZKurtosis = mean(combinedZKurtosis);
    stdCombinedZKurtosis = std(combinedZKurtosis);
    maxCombinedZKurtosis = max(combinedZKurtosis);
    minCombinedZKurtosis = min(combinedZKurtosis);
    
    fprintf(summaryFile, 'RMS Probe Percent Difference:\n Mean:\t%3.2f%%\n Std Dev:%3.2f%%\t\n Max:\t%3.2f%%\n Min:\t%3.2f%%\n',...
        meanProbeRMSPercentDiff, stdProbeRMSPercentDiff, maxProbeRMSPercentDiff, minProbeRMSPercentDiff);
    fprintf(summaryFile, 'X Skewness Probe:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanProbeXSkewness, stdProbeXSkewness, maxProbeXSkewness, minProbeXSkewness);
    fprintf(summaryFile, 'Y Skewness Probe:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanProbeYSkewness, stdProbeYSkewness, maxProbeYSkewness, minProbeYSkewness);
    fprintf(summaryFile, 'Z Skewness Probe:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanProbeZSkewness, stdProbeZSkewness, maxProbeZSkewness, minProbeZSkewness);
    fprintf(summaryFile, 'X Kurtosis Probe:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanProbeXKurtosis, stdProbeXKurtosis, maxProbeXKurtosis, minProbeXKurtosis);
    fprintf(summaryFile, 'Y Kurtosis Probe:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanProbeYKurtosis, stdProbeYKurtosis, maxProbeYKurtosis, minProbeYKurtosis);
    fprintf(summaryFile, 'Z Kurtosis Probe:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanProbeZKurtosis, stdProbeZKurtosis, maxProbeZKurtosis, minProbeZKurtosis);
    
    fprintf(summaryFile, 'RMS Reference Percent Difference:\n Mean:\t%3.2f%%\n Std Dev:%3.2f%%\t\n Max:\t%3.2f%%\n Min:\t%3.2f%%\n',...
        meanRefRMSPercentDiff, stdRefRMSPercentDiff, maxRefRMSPercentDiff, minRefRMSPercentDiff);
    fprintf(summaryFile, 'X Skewness Ref:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanRefXSkewness, stdRefXSkewness, maxRefXSkewness, minRefXSkewness);
    fprintf(summaryFile, 'Y Skewness Ref:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanRefYSkewness, stdRefYSkewness, maxRefYSkewness, minRefYSkewness);
    fprintf(summaryFile, 'Z Skewness Ref:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanRefZSkewness, stdRefZSkewness, maxRefZSkewness, minRefZSkewness);
    fprintf(summaryFile, 'X Kurtosis Ref:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanRefXKurtosis, stdRefXKurtosis, maxRefXKurtosis, minRefXKurtosis);
    fprintf(summaryFile, 'Y Kurtosis Ref:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanRefYKurtosis, stdRefYKurtosis, maxRefYKurtosis, minRefYKurtosis);
    fprintf(summaryFile, 'Z Kurtosis Ref:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanRefZKurtosis, stdRefZKurtosis, maxRefZKurtosis, minRefZKurtosis);
    
    fprintf(summaryFile, 'RMS Combined Percent Difference:\n Mean:\t%3.2f%%\n Std Dev:%3.2f%%\t\n Max:\t%3.2f%%\n Min:\t%3.2f%%\n',...
        meanCombinedRMSPercentDiff, stdCombinedRMSPercentDiff, maxCombinedRMSPercentDiff, minCombinedRMSPercentDiff);
    fprintf(summaryFile, 'X Skewness Combined:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanCombinedXSkewness, stdCombinedXSkewness, maxCombinedXSkewness, minCombinedXSkewness);
    fprintf(summaryFile, 'Y Skewness Combined:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanCombinedYSkewness, stdCombinedYSkewness, maxCombinedYSkewness, minCombinedYSkewness);
    fprintf(summaryFile, 'Z Skewness Combined:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanCombinedZSkewness, stdCombinedZSkewness, maxCombinedZSkewness, minCombinedZSkewness);
    fprintf(summaryFile, 'X Kurtosis Combined:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanCombinedXKurtosis, stdCombinedXKurtosis, maxCombinedXKurtosis, minCombinedXKurtosis);
    fprintf(summaryFile, 'Y Kurtosis Combined:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanCombinedYKurtosis, stdCombinedYKurtosis, maxCombinedYKurtosis, minCombinedYKurtosis);
    fprintf(summaryFile, 'Z Kurtosis Combined:\n Mean:\t%3.2f\n Std Dev:%3.2f\t\n Max:\t%3.2f\n Min:\t%3.2f\n',...
        meanCombinedZKurtosis, stdCombinedZKurtosis, maxCombinedZKurtosis, minCombinedZKurtosis);

    maxValue = max([probePredictedRMS; probeMeasuredRMS]);
    minValue = min([probePredictedRMS; probeMeasuredRMS]);
    
    figure(4);
    plot(probePredictedRMS, probeMeasuredRMS, 'k.');
    hold on;
    plot([minValue maxValue], [minValue, maxValue], 'k');
    xlabel('Predicted TRE RMS (mm)', 'fontsize',14);
    ylabel('Simulated TRE RMS (mm)', 'fontsize',14);
    set(gca,'fontsize',14);
    title('Comparison of Predicted and Simulated TRE RMS for Probe', 'fontsize',14);
    hold off;
    axis([minValue, maxValue, minValue, maxValue]);
    figurefilename = sprintf('PredvMeas_%s_Probe', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);

    maxValue = max([refPredictedRMS; refMeasuredRMS]);
    minValue = min([refPredictedRMS; refMeasuredRMS]);

    figure(5);
    plot(refPredictedRMS, refMeasuredRMS, 'k.');
    hold on;
    plot([minValue maxValue], [minValue, maxValue], 'k');
    xlabel('Predicted TRE RMS (mm)', 'fontsize',14);
    ylabel('Simulated TRE RMS (mm)', 'fontsize',14);
    set(gca,'fontsize',14);
    title('Comparison of Predicted and Simulated TRE RMS for Reference', 'fontsize',14);
    hold off;
    axis([minValue, maxValue, minValue, maxValue]);
    figurefilename = sprintf('PredvMeas_%s_Reference', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);
    
    maxValue = max([combinedPredictedRMS; combinedMeasuredRMS]);
    minValue = min([combinedPredictedRMS; combinedMeasuredRMS]);

    figure(6);
    plot(combinedPredictedRMS, combinedMeasuredRMS, 'k.');
    hold on;
    plot([minValue maxValue], [minValue, maxValue], 'k');
    xlabel('Predicted TRE RMS (mm)', 'fontsize',14);
    ylabel('Simulated TRE RMS (mm)', 'fontsize',14);
    set(gca,'fontsize',14);
    title('Comparison of Predicted and Simulated TRE RMS Combined', 'fontsize',14);
    hold off;
    axis([minValue, maxValue, minValue, maxValue]);
    figurefilename = sprintf('PredvMeas_%s_Reference', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);

    
    fclose(outFile);
    fclose(summaryFile);

    cd ..
end

fclose('all');
cd 'E:\docs\research\phd\experiments\TRE Experiments'
return;
%%
combined.Actual.stats.Sigma
combined.Meas.stats.Sigma

figure(7);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', combined.Meas.stats.Sigma, 'b', 1);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', combined.Actual.stats.Sigma, 'g', 1);

figure(8);
plotStochastic( combined.Meas.error', [0 0 0]', combined.Meas.stats.Sigma, [0 0 0]', combined.Actual.stats.Sigma );

figure(9);
plot3(probe.Actual.mrk(:,1), probe.Actual.mrk(:,2), probe.Actual.mrk(:,3),'b.');
hold on;
plot3(probe.Actual.tip(:,1), probe.Actual.tip(:,2), probe.Actual.tip(:,3),'b+');
plot3(ref.Actual.mrk(:,1), ref.Actual.mrk(:,2), ref.Actual.mrk(:,3),'k.');
plot3(ref.Actual.tip(:,1), ref.Actual.tip(:,2), ref.Actual.tip(:,3),'ko');
hold off;

figure(10);
subplot(2,2,1);
plot3(probe.Meas.tip(:,1),probe.Meas.tip(:,2),probe.Meas.tip(:,3), '.b');
hold on;
plot3(ref.Meas.tip(:,1),ref.Meas.tip(:,2),ref.Meas.tip(:,3), '.k');
test = (refXfrm.R*combined.Meas.tip')' + repmat(refXfrm.pos, N,1);
%plot3(combined.Meas.tip(:,1),combined.Meas.tip(:,2),combined.Meas.tip(:,3), '.r');
plot3(test(:,1),test(:,2),test(:,3), '.r');
hold off;
subplot(2,2,2);
plot3(probe.Meas.tip(:,1),probe.Meas.tip(:,2),probe.Meas.tip(:,3), '.b');
subplot(2,2,3);
plot3(ref.Meas.tip(:,1),ref.Meas.tip(:,2),ref.Meas.tip(:,3), '.k');
subplot(2,2,4);
plot3(test(:,1),test(:,2),test(:,3), '.r');

figure(11);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', combined.Meas.stats.Sigma, 'r', 1);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', probe.Meas.stats.Sigma, 'b', 1);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', ref.Meas.stats.Sigma', 'k', 1);

figure(12);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', combined.Actual.stats.Sigma, 'r', 1);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', probe.Actual.stats.Sigma, 'b', 1);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', ref.Actual.stats.Sigma, 'k', 1);

figure(13);
test = [probe.Meas.error ref.Meas.error];
testCov6 = cov(test);
testCovA = testCov6(1:3, 1:3) + testCov6(4:6, 4:6) 
testCovB = testCov6(1:3, 1:3) + testCov6(4:6, 4:6) + 2* testCov6(1:3,4:6)

refXfrm.R'*combined.Meas.stats.Sigma*refXfrm.R
refXfrm.R'*combined.Actual.stats.Sigma*refXfrm.R

figure(14);
plotStochastic( ref.Meas.error', [0 0 0]', ref.Meas.stats.Sigma, [0 0 0]', ref.Actual.stats.Sigma );