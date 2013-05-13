%% simulate a bunch of error xfrms.
close all;
clear all;
%email_setup;
%% set up the monte carlo parameters.
N = 100000;  %no. of samples for each monte carlo.
M = 1000;     %no. of sims to run.
sensor = '6D';
scale = 4:4:20;

% on iron.irus.robarts.ca:
%cd E:\awiles\data\TTEMonteCarlo
% on home Mac:
cd E:\data\TTEMonteCarlo
%cd E:\docs\research\phd\experiments\TTEMonteCarlo
switch(sensor)
    case '5D'
        load target5D.mat;
        ftefilename = '5DCov-MonteCarlo.csv';
        xfrmfilename = '5DXfrms-MonteCarlo.csv';
        name =  'TTENeedle';
    case '6D'
        load target.mat;
        ftefilename = '6DCov-MonteCarlo.csv';
        xfrmfilename = '6DXfrms-MonteCarlo.csv';
        name = 'TTERadolfzell'
    otherwise
        error('Invalid sensor type given: %s', sensor);
end

for targetid = 1:3
    for sampleNo=5
        for nScale = 1:length(scale)
            %% make sure we are in the correct directory.
            % on iron.irus.robarts.ca:
            %cd E:\awiles\data\TTEMonteCarlo
            % on home Mac:
            cd E:\data\TTEMonteCarlo
            %cd E:\docs\research\phd\experiments\TTEMonteCarlo

            %% set the casename
            casename = sprintf('%s-tar%02d-xfrm%02d-r2scale%02d', name, targetid, sampleNo, scale(nScale));

            %% get the FTE matrix.
            fte = readFTECovariance(ftefilename,sampleNo);
            %fte = scale(nScale) * fte;

            %% get the transforms at which we are working.
            truexfrm = readXfrmFile(xfrmfilename, sampleNo);

            %% get the random targets:
            r0 = r(targetid,:);
            r0 = scale(nScale) * r0;

            %% run the simulation.
            mkdir(casename);
            cd(casename);
            sumdata = zeros(M,9); %rms % diff, wishart local cov, wishart global cov, wishart local meanandcov, wishart global meanandcov

            starttime = cputime;
            for i = 1:M
                currentcasename = sprintf('%s-%04d', casename, i);
                [data, stats] = procTTEMonteCarlo(truexfrm, fte, r0, N, [currentcasename '.txt']);
                sumdata(i, :) = [stats.rmspercentdiff, stats.wishart.local.covariance,...
                    stats.wishart.global.covariance, stats.wishart.local.meanandcov,...
                    stats.wishart.global.meanandcov,stats.wishartFirstOrder.local.covariance,...
                    stats.wishartFirstOrder.global.covariance, stats.wishartFirstOrder.local.meanandcov,...
                    stats.wishartFirstOrder.global.meanandcov];
                fprintf('%s % 4d -- RMS %% Diff: % 3.2f +- % 2.2f   \n',casename, i, mean(sumdata(1:i,1)), std(sumdata(1:i,1)));
                fprintf('Wishart Higher Order %%: % 2.2f%%, % 2.2f%%, % 2.2f%%, % 2.2f%%\n', 100*sum(sumdata(:,2))/i, ...
                    100*sum(sumdata(:,3))/i, 100*sum(sumdata(:,4))/i, 100*sum(sumdata(:,5))/i);
                fprintf('Wishart 1st Order    %%: % 2.2f%%, % 2.2f%%, % 2.2f%%, % 2.2f%%\n\n', 100*sum(sumdata(:,6))/i, ...
                    100*sum(sumdata(:,7))/i, 100*sum(sumdata(:,8))/i, 100*sum(sumdata(:,9))/i);
                save(currentcasename, 'data', 'stats');

                % compute estimated completion time.
                avgTime = (cputime - starttime)/i; %in seconds.
                estTimeToFinish = avgTime*(M-i); %in seconds.
                estHours = floor(estTimeToFinish/3600);
                estMinutes = floor((estTimeToFinish - 3600*estHours)/60);

                fprintf('Completed %3.1f %%... Avg. time per run: %4.2f seconds. Finish in approx. %d hours %d min \n',...
                    (100*i/M), avgTime, estHours, estMinutes);
            end
            save(casename, 'sumdata');
            flog = fopen([casename '.txt'], 'wt');
            fprintf(flog, 'No. of Runs %d\n -- RMS %% Diff: % 3.2f +- % 2.2f\n', i, mean(sumdata(:,1)), std(sumdata(:,1)));
            fprintf(flog,'Wishart Higher Order %%: % 2.2f%%, % 2.2f%%, % 2.2f%%, % 2.2f%%\n', 100*sum(sumdata(:,2))/i, ...
                100*sum(sumdata(:,3))/i, 100*sum(sumdata(:,4))/i, 100*sum(sumdata(:,5))/i);
            fprintf(flog,'Wishart 1st Order    %%: % 2.2f%%, % 2.2f%%, % 2.2f%%, % 2.2f%%\n', 100*sum(sumdata(:,6))/i, ...
                100*sum(sumdata(:,7))/i, 100*sum(sumdata(:,8))/i, 100*sum(sumdata(:,9))/i);
            fclose(flog);
            % sendmail('adwiles@gmail.com', ['TTE Simulation Complete: ', casename], 'Report is attached.', [casename '.txt']);
            cd ..
            %% plot the data.
            %plotTTEMonteCarlo(data, stats);

            %% create a 3D sensor error plot.

            % note FTE from above.
            %step = 10;
            %nEllipsoids = 10;

            %figure(gcf+1);
            %plotErrorPropagation(fte, step, nEllipsoids, truexfrm);

            %figure(gcf+1);

            %plotAxes(nEllipsoids*step + step/2, [2,3,1], truexf
        end
    end
end
