function [data, stats] = procTTEMonteCarlo(truexfrm, fte, r0, N, logfilename)
% procTTEMonteCarlo
% truexfrm -- transform of the sensor associated with the FTE.
% fte -- fiducial tracking error covariance matrix.
% r0 -- target location in local sensor coordinate frame.
% N -- number of samples in the Monte Carlo simulation.

%% find the true global position of the target.
rg = getXfrmPointQuat(truexfrm,r0); %get the target in true global space.

%% compute the stats from the TTE formula.
[stats.theory.local.mu, stats.theory.local.cov, stats.theory.local.rms] = calcTTE(fte,r0);
[stats.theoryFirstOrder.local.mu, stats.theoryFirstOrder.local.cov, stats.theoryFirstOrder.local.rms] = calcTTE(fte,r0,1);
% rotate the stats into the global frame.
[stats.theory.global.mu, stats.theory.global.cov] = rotateCovQuat(stats.theory.local.mu,...
    stats.theory.local.cov, truexfrm.rot);
stats.theory.global.rms = stats.theory.local.rms;
[stats.theoryFirstOrder.global.mu, stats.theoryFirstOrder.global.cov] = rotateCovQuat(stats.theoryFirstOrder.local.mu,...
    stats.theoryFirstOrder.local.cov, truexfrm.rot);
stats.theoryFirstOrder.global.rms = stats.theoryFirstOrder.local.rms;

%% initialize variables and run the simulations.
rlocal = zeros(N,3);
ttelocal = zeros(N,3);
xfrmstack = zeros(N,7);

r = zeros(N,3);
tte = zeros(N,3);

for i = 1:N
    xfrm = getRandomXfrm(fte);
    % store the transform.
    xfrmstack(i,1:3) = xfrm.pos;
    xfrmstack(i,4:6) = xfrm.rot(2:4);
    xfrmstack(i,7) = xfrm.rot(1);
    % compute the transform error.
    rlocal(i,:) = getXfrmPointQuat(xfrm,r0);
    r(i,:) = getXfrmPointQuat(truexfrm,rlocal(i,:));
    ttelocal(i,:) = rlocal(i,:) - r0;
    tte(i,:) = r(i,:) - rg;
end


xfrmcov = cov(xfrmstack);
stats.xfrmcov = xfrmcov;
data.xfrmstack = xfrmstack;
data.tte = tte;
data.ttelocal = ttelocal;
%% compute the stats.
stats.simulated.global.mu = mean(tte);
stats.simulated.global.cov = cov(tte);
stats.simulated.global.rms = sqrt(mean(sum(tte.^2,2)));

stats.simulated.local.mu = mean(ttelocal);
stats.simulated.local.cov = cov(ttelocal);
stats.simulated.local.rms = sqrt(mean(sum(ttelocal.^2,2)));

%% perform the Wishart Hypothesis tests - higher order model.
[chi2stat, chi2prob, dof] = WishartHypTest(stats.theory.local.mu,...
    stats.theory.local.cov, ttelocal, 'covariance');
stats.wishart.local.covariance = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(stats.theory.global.mu,...
    stats.theory.global.cov, tte, 'covariance');
stats.wishart.global.covariance = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(stats.theory.local.mu,...
    stats.theory.local.cov, ttelocal, 'meanandcov');
stats.wishart.local.meanandcov = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(stats.theory.global.mu,...
    stats.theory.global.cov, tte, 'meanandcov');
stats.wishart.global.meanandcov = checkChi2Test(chi2stat, dof);

%% perform the Wishart Hypothesis tests - first order model.
[chi2stat, chi2prob, dof] = WishartHypTest(stats.theoryFirstOrder.local.mu,...
    stats.theoryFirstOrder.local.cov, ttelocal, 'covariance');
stats.wishartFirstOrder.local.covariance = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(stats.theoryFirstOrder.global.mu,...
    stats.theoryFirstOrder.global.cov, tte, 'covariance');
stats.wishartFirstOrder.global.covariance = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(stats.theoryFirstOrder.local.mu,...
    stats.theoryFirstOrder.local.cov, ttelocal, 'meanandcov');
stats.wishartFirstOrder.local.meanandcov = checkChi2Test(chi2stat, dof);

[chi2stat, chi2prob, dof] = WishartHypTest(stats.theoryFirstOrder.global.mu,...
    stats.theoryFirstOrder.global.cov, tte, 'meanandcov');
stats.wishartFirstOrder.global.meanandcov = checkChi2Test(chi2stat, dof);


%% get the percent difference for the RMS.
stats.rmspercentdiff = 100*(stats.simulated.global.rms ...
    - stats.theory.global.rms)/stats.theory.global.rms;
%% get the percent difference for the RMS.
stats.rmspercentdiffFirstOrder = 100*(stats.simulated.global.rms ...
    - stats.theoryFirstOrder.global.rms)/stats.theoryFirstOrder.global.rms;

%% write stats to log.
if(nargin > 4)
    flog = fopen(logfilename, 'wt');

    fprintf(flog, 'TTE Simulation Log\n\n');
    fprintf(flog, 'Sensor Transformation\n\n');
    fprintf(flog, 'Pos:  [ % g  % g  % g ]\n', truexfrm.pos);
    fprintf(flog, 'Quat: [ % g  % g  % g  % g ]\n\n', truexfrm.rot);
    fprintf(flog, 'Local Target:  [ % g  % g  % g ]\n', r0);
    fprintf(flog, 'Global Target:  [ % g  % g  % g ]\n\n', rg);
    fprintf(flog, 'FTE Input Covariance Matrix: \n');
    if(size(fte,1) == 5)
        for i = 1:5
            fprintf(flog, '|% g\t% g\t% g\t% g\t% g|\n', fte(i,:));
        end
    elseif(size(fte,1) == 6)
        for i = 1:6
            fprintf(flog, '|% g\t% g\t% g\t% g\t% g\t% g|\n', fte(i,:));
        end
    else
        error('procTTEMonteCarlo::Someweird size of the FTE.');
    end
    fprintf(flog, '\nFTE Actual Covariance Matrix (with q0 as the 7th variable): \n');
    for i = 1:7
        fprintf(flog, '|% g\t% g\t% g\t% g\t% g\t% g\t% g|\n', xfrmcov(i,:));
    end
    fprintf(flog, '\nStatistics -- Higher Order\n\nMean -- Local Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\t%% Diff\n');
    if(size(fte,1) == 5)
        start =3;
    elseif(size(fte,1) == 6)
        start = 1;
    else
        error('procTTEMonteCarlo::Someweird size of the FTE.');
    end
    for i = start:3
        fprintf(flog, '% g\t% g\t% 2.2f%%\n', stats.theory.local.mu(i), stats.simulated.local.mu(i),...
            100*((stats.simulated.local.mu(i)-stats.theory.local.mu(i))/stats.theory.local.mu(i)));
    end
    fprintf(flog, '\nMean -- Global Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\t%% Diff\n');
    for i = 1:3
        fprintf(flog, '% g\t% g\t% 2.2f%%\n', stats.theory.global.mu(i), stats.simulated.global.mu(i),...
            100*((stats.simulated.global.mu(i)-stats.theory.global.mu(i))/stats.theory.global.mu(i)));
    end

    percdiff1 = 100*(stats.simulated.local.cov - stats.theory.local.cov)./stats.theory.local.cov;
    percdiff2 = 100*(stats.simulated.global.cov - stats.theory.global.cov)./stats.theory.global.cov;
    fprintf(flog, '\nCovariance -- Local Frame\n');
    fprintf(flog, 'Theory\t\t\t\t\t\tSimulation\t\t\t\t\t\t%% Diff\n');
    for i = 1:3
        fprintf(flog, '|% g\t% g\t% g |  |% g\t% g\t% g |  |% 2.2f%%\t% 2.2f%%\t% 2.2f%% |\n',...
            stats.theory.local.cov(i,:), stats.simulated.local.cov(i,:), percdiff1(i,:));
    end

    fprintf(flog, '\nCovariance -- Global Frame\n');
    fprintf(flog, 'Theory\t\t\t\t\t\tSimulation\t\t\t\t\t\t%% Diff\n');
    for i = 1:3
        fprintf(flog, '|% g\t% g\t% g |  |% g\t% g\t% g |  |% 2.2f%%\t% 2.2f%%\t% 2.2f%% |\n',...
            stats.theory.global.cov(i,:), stats.simulated.global.cov(i,:), percdiff2(i,:));
    end

    fprintf(flog, '\nRMS -- Local Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\n');
    fprintf(flog, '% g\t% g\n', stats.theory.local.rms, stats.simulated.local.rms);

    fprintf(flog, '\nRMS -- Global Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\t%% Diff\n');
    fprintf(flog, '% g\t% g\t% 2.2f%%\n', stats.theory.global.rms, stats.simulated.global.rms, stats.rmspercentdiff);

    fprintf(flog, '\nWishart Hypothesis Tests\n');
    fprintf(flog, '\tCovariance\tMean and Cov.\n');
    fprintf(flog, 'Local\t%d\t\t%d\n', stats.wishart.local.covariance, stats.wishart.local.meanandcov);
    fprintf(flog, 'Global\t%d\t\t%d\n\n\n', stats.wishart.global.covariance, stats.wishart.global.meanandcov);
    
    fprintf(flog, '\nStatistics -- First Order \n\nMean -- Local Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\t%% Diff\n');
    if(size(fte,1) == 5)
        start =3;
    elseif(size(fte,1) == 6)
        start = 1;
    else
        error('procTTEMonteCarlo::Someweird size of the FTE.');
    end
    for i = start:3
        fprintf(flog, '% g\t% g\t% 2.2f%%\n', stats.theoryFirstOrder.local.mu(i), stats.simulated.local.mu(i),...
            100*((stats.simulated.local.mu(i)-stats.theoryFirstOrder.local.mu(i))/stats.theoryFirstOrder.local.mu(i)));
    end
    fprintf(flog, '\nMean -- Global Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\t%% Diff\n');
    for i = 1:3
        fprintf(flog, '% g\t% g\t% 2.2f%%\n', stats.theoryFirstOrder.global.mu(i), stats.simulated.global.mu(i),...
            100*((stats.simulated.global.mu(i)-stats.theoryFirstOrder.global.mu(i))/stats.theoryFirstOrder.global.mu(i)));
    end

    percdiff1 = 100*(stats.simulated.local.cov - stats.theoryFirstOrder.local.cov)./stats.theoryFirstOrder.local.cov;
    percdiff2 = 100*(stats.simulated.global.cov - stats.theoryFirstOrder.global.cov)./stats.theoryFirstOrder.global.cov;
    fprintf(flog, '\nCovariance -- Local Frame\n');
    fprintf(flog, 'Theory\t\t\t\t\t\tSimulation\t\t\t\t\t\t%% Diff\n');
    for i = 1:3
        fprintf(flog, '|% g\t% g\t% g |  |% g\t% g\t% g |  |% 2.2f%%\t% 2.2f%%\t% 2.2f%% |\n',...
            stats.theoryFirstOrder.local.cov(i,:), stats.simulated.local.cov(i,:), percdiff1(i,:));
    end

    fprintf(flog, '\nCovariance -- Global Frame\n');
    fprintf(flog, 'Theory\t\t\t\t\t\tSimulation\t\t\t\t\t\t%% Diff\n');
    for i = 1:3
        fprintf(flog, '|% g\t% g\t% g |  |% g\t% g\t% g |  |% 2.2f%%\t% 2.2f%%\t% 2.2f%% |\n',...
            stats.theoryFirstOrder.global.cov(i,:), stats.simulated.global.cov(i,:), percdiff2(i,:));
    end

    fprintf(flog, '\nRMS -- Local Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\n');
    fprintf(flog, '% g\t% g\n', stats.theoryFirstOrder.local.rms, stats.simulated.local.rms);

    fprintf(flog, '\nRMS -- Global Frame\n');
    fprintf(flog, 'Theory\t\tSimulation\t%% Diff\n');
    fprintf(flog, '% g\t% g\t% 2.2f%%\n', stats.theoryFirstOrder.global.rms, stats.simulated.global.rms, stats.rmspercentdiffFirstOrder);

    fprintf(flog, '\nWishart Hypothesis Tests\n');
    fprintf(flog, '\tCovariance\tMean and Cov.\n');
    fprintf(flog, 'Local\t%d\t\t%d\n', stats.wishartFirstOrder.local.covariance, stats.wishartFirstOrder.local.meanandcov);
    fprintf(flog, 'Global\t%d\t\t%d\n', stats.wishartFirstOrder.global.covariance, stats.wishartFirstOrder.global.meanandcov);

    fclose(flog);
end