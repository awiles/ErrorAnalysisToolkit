function testcase = computeStats(testcase)

%get summary statistics.
testcase.Meas.stats.mu = mean(testcase.Meas.error);
testcase.Meas.stats.Sigma = cov(testcase.Meas.error);
testcase.Meas.stats.RMS = sqrt(mean(sum(testcase.Meas.error.^2,2)));
testcase.Meas.stats.max = max(testcase.Meas.error);
testcase.Meas.stats.min = min(testcase.Meas.error);
testcase.Meas.stats.absmax = max(abs(testcase.Meas.error));
testcase.Meas.stats.skewness.x = getSkewness(testcase.Meas.tip(:,1));
testcase.Meas.stats.skewness.y = getSkewness(testcase.Meas.tip(:,2));
testcase.Meas.stats.skewness.z = getSkewness(testcase.Meas.tip(:,3));
testcase.Meas.stats.kurtosis.x = getKurtosis(testcase.Meas.tip(:,1));
testcase.Meas.stats.kurtosis.y = getKurtosis(testcase.Meas.tip(:,2));
testcase.Meas.stats.kurtosis.z = getKurtosis(testcase.Meas.tip(:,3));