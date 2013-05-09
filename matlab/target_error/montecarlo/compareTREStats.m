function testResult = compareTREStats(probe, ref, combined)
%**************************************************************************
% Compare the TRE from Theory and Monte Carlo Simulation.
%
%       Written by Andrew Wiles, May 23, 2007
%
%**************************************************************************

%% compute and check the probe stats.
%************RMS Percent Diff********************%
testResult.probe.RMS.Meas = probe.Meas.stats.RMS;
testResult.probe.RMS.Theory = probe.Actual.stats.RMS;
testResult.probe.RMS.Diff =  probe.Actual.stats.RMS - probe.Meas.stats.RMS;
testResult.probe.RMS.PercentDiff =...
    100*(testResult.probe.RMS.Diff/probe.Meas.stats.RMS);

if(~probe.bHomogenous)
    testResult.probeAvgFLE.RMS.Meas = probe.Meas.stats.RMS;
    testResult.probeAvgFLE.RMS.Theory = probe.AvgFLE.stats.RMS;
    testResult.probeAvgFLE.RMS.Diff = probe.AvgFLE.stats.RMS - probe.Meas.stats.RMS;
    testResult.probeAvgFLE.RMS.PercentDiff =...
        100*(testResult.probeAvgFLE.RMS.Diff/probe.Meas.stats.RMS);
end
%************Wishart Test************************%
[chi2stat, chi2prob, dof] = WishartHypTest(probe.Actual.stats.mu,...
    probe.Actual.stats.Sigma, probe.Meas.error, 'covariance');
testResult.probe.covariance = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(probe.Actual.stats.mu,...
    probe.Actual.stats.Sigma, probe.Meas.error, 'meanandcov');
probe.Actual.stats.chi2stat = chi2stat;
testResult.probe.meanandcov = checkChi2Test(chi2stat, dof);

if(~probe.bHomogenous)
    [chi2stat, chi2prob, dof] = WishartHypTest(probe.AvgFLE.stats.mu,...
        probe.AvgFLE.stats.Sigma, probe.Meas.error, 'covariance');
    testResult.probeAvgFLE.covariance = checkChi2Test(chi2stat, dof);

    [chi2stat, chi2prob, dof] = WishartHypTest(probe.AvgFLE.stats.mu,...
        probe.AvgFLE.stats.Sigma, probe.Meas.error, 'meanandcov');
    probe.AvgFLE.stats.chi2stat = chi2stat;
    testResult.probeAvgFLE.meanandcov = checkChi2Test(chi2stat, dof);
end

if(nargin > 1) % if reference test is being run.
    %% TODO: update for non-homogeonous results.
    %% check the reference stats.
    %************RMS Percent Diff********************%
    testResult.ref.RMS.Meas = ref.Meas.stats.RMS;
    testResult.ref.RMS.Theory = ref.Actual.stats.RMS;
    testResult.ref.RMS.Diff = ref.Meas.stats.RMS - ref.Actual.stats.RMS;
    testResult.ref.RMS.PercentDiff =...
        100*(testResult.ref.RMS.Diff/ref.Meas.stats.RMS);

    %************Wishart Test************************%
    [chi2stat, chi2prob, dof] = WishartHypTest(ref.Actual.stats.mu,...
        ref.Actual.stats.Sigma, ref.Meas.error, 'covariance');
    testResult.ref.covariance = checkChi2Test(chi2stat, dof);

    [chi2stat, chi2prob, dof] = WishartHypTest(ref.Actual.stats.mu,...
        ref.Actual.stats.Sigma, ref.Meas.error, 'meanandcov');
    ref.Actual.stats.chi2stat = chi2stat;
    testResult.ref.meanandcov = checkChi2Test(chi2stat, dof);

    %% check the combined stats.
    %************RMS Percent Diff********************%
    testResult.combined.RMS.Meas = combined.Meas.stats.RMS;
    testResult.combined.RMS.Theory = combined.Actual.stats.RMS;
    testResult.combined.RMS.Diff = combined.Meas.stats.RMS - combined.Actual.stats.RMS;
    testResult.combined.RMS.PercentDiff =...
        100*(testResult.combined.RMS.Diff/combined.Meas.stats.RMS);

    %************Wishart Test************************%
    [chi2stat, chi2prob, dof] = WishartHypTest(combined.Actual.stats.mu,...
        combined.Actual.stats.Sigma, combined.Meas.error, 'covariance');
    testResult.combined.covariance = checkChi2Test(chi2stat, dof);

    [chi2stat, chi2prob, dof] = WishartHypTest(combined.Actual.stats.mu,...
        combined.Actual.stats.Sigma, combined.Meas.error, 'meanandcov');
    combined.Actual.stats.chi2stat = chi2stat;
    testResult.combined.meanandcov = checkChi2Test(chi2stat, dof);
end

%%
function result = checkChi2Test(chi2stat, dof)
% assume alpha = 0.05.
if( chi2stat < 0)
    error('Invalid chi2 statistic returned.');
end

switch(dof)
    case 6
        if( chi2stat < 12.59)
            result = 1;
        else
            result = 0;
        end
    case 9
        if( chi2stat < 16.92)
            result = 1;
        else
            result = 0;
        end
    otherwise
        error('Invalid dof returned.  Check the data.');
end

