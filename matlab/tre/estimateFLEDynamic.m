%function estimateFLEDynamic(casename, casestudy, winsize)
function [data, stats] = estimateFLEDynamic(tooldesign, flemodel, toolpath, winsize, testnum, bPlot)
%% simulator for estimating FLE.
%cd E:\docs\research\phd\experiments\FLEPrediction\simulated-tool-paths
casename = sprintf('%s-%s-%s-%03d-%04d', tooldesign.name, flemodel.name,...
    toolpath.name, winsize, testnum);
outfilename = sprintf('%s.txt', casename);
outfilename2 = sprintf('%sLaTeX.txt',casename);
fid2 = fopen(outfilename2, 'wt');
fid = fopen(outfilename, 'wt');
fprintf(fid, 'Casename: %s\n\n', casename);
fprintf(fid2, 'Casename: %s\n\n', casename);

refmrk = tooldesign.refmrk;
truexfrm = toolpath.xfrm;

nMrks = size(refmrk,1);
N = length(truexfrm); % xfrm is from the path-%02d.mat file.

% some old stuff that we'll keep here for now.
covmethod = 1;
pamethod = 1; %method used to rotate into principal axes: 1: xfrmed positions, 2: true positions.

%% fle model parameters.
%model = 'radial';
%rate = [5.7e-9, 5.7e-9, 5.1e-8];
%% characteristic matrix.
% the following call was replaced with one outside of the simulation.
%[U, L, V0, x, A, Ainv] = xfrmToPA(refmrk);
A = tooldesign.A;
Ainv = tooldesign.Ainv;
V0 = tooldesign.V0;

%% set up the Monte Carlo simulation.
frewin = cell(nMrks);
frecovwin = cell(nMrks);
fleRMS = cell(nMrks);
fleValues = cell(nMrks);
fleEstValues = cell(nMrks);
fleEigValues = cell(nMrks);
fleEstEigValues = cell(nMrks);
fremean = cell(nMrks);
frecov = cell(nMrks);
frecovmean = cell(nMrks);
frenoisycov = cell(nMrks);
fresmoothcov = cell(nMrks);
RPA = cell(nMrks);
for i = 1:nMrks
    % perform some preallocations.
    frewin{i} = zeros(winsize,3);
    frecovwin{i} = zeros(winsize,6);
    fleRMS{i} = zeros(N,1);
    fleValues{i} = zeros(6,N);
    fleEstValues{i} = zeros(6,N);
    fleEigValues{i} = zeros(3,N);
    fleEstEigValues{i} = zeros(3,N);
    fremean{i} = [0 0 0];
    frecov{i} = zeros(3);
    frecovmean{i} = [0 0 0 0 0 0]; %for the moving average.
    frenoisycov{i} = zeros(N,6);
    fresmoothcov{i} = zeros(N,6);
end
measfreRMS = zeros(N, nMrks);

fleestRMS = zeros(N,nMrks);
meanFLERMS = zeros(N,1);
meanFLEEigValues = zeros(3,N);
%RRMS = zeros(N,nMrks);
frecovRMS = zeros(N,nMrks);
fresmoothRMS = zeros(N,nMrks);
%% run the loop.
starttime = cputime;
for i = 1:N
    % get the marker positions.
    mrk = (truexfrm{i}.R * refmrk')' + repmat(truexfrm{i}.pos, size(refmrk,1), 1);

    %compute centroid position.
    cent = mean(mrk);
    %set up the fle.
    if(flemodel.heteroscedastic)
        fle = getQuadraticFLE(mrk, flemodel.rate, flemodel.model);
    else
        fle = getQuadraticFLE(cent, flemodel.rate, flemodel.model);
    end
    %add the noise.
    measmrk = mrk + getMeasNoise(fle, nMrks);

    %*****************************************
    %compute the transform and xfrm positions.
    %*****************************************
    % using SVD:
    xfrm = getRigidXfrmSVD( refmrk, measmrk );
    xfrmmrk = (xfrm.rot * refmrk')' + repmat(xfrm.pos, nMrks, 1);

    % using Horn:
    %xfrm = getRigidXfrm( refmrk, measmrk );
    %xfrmmrk = (quat2rm(xfrm.rot) * refmrk')' + repmat(xfrm.pos, nMrks, 1);

    switch pamethod
        case 1
            %[U, L, V, x, A, Ainv] = xfrmToPA(xfrmmrk);
            %[U, L, V] = xfrmToPA(xfrmmrk);
            %[V, D] = getPA(xfrmmrk);
            V = xfrm.rot*V0;
            %V = V0*xfrm.rot';
        case 2
            %[U, L, V, x, A, Ainv] = xfrmToPA(mrk);
            %[U, L, V] = xfrmToPA(mrk);
            [V, D] = getPA(mrk);
        otherwise
            error('Invalid method to rotate into PA given.');
    end
    %compute the FRE.
    fre = xfrmmrk - measmrk;
    meanFLE = zeros(3);
    for j=1:nMrks
        switch covmethod
            case 1
                %update the mean.
                fremean{j} = fremean{j} + 1/winsize* (fre(j,:) - frewin{j}(end,:));
                %update the covariance
                frecov{j} = frecov{j} + 1/(winsize-1) * (fre(j,:)'*fre(j,:) - frewin{j}(end,:)'*frewin{j}(end,:));
                %update the window.
                frewin{j} = circshift(frewin{j},1);
                frewin{j}(1,:) = fre(j,:);
            case 2
                %update the window.
                frewin{j} = circshift(frewin{j},1);
                frewin{j}(1,:) = fre(j,:);
                fremean{j} = mean(frewin{j});
                frecov{j} = cov(frewin{j}); %TODO: replace this with more efficient update code.
            otherwise
                error('Invalid covariance method given.');
        end

        %rotate R into PA.
        %RPA{j} = V' * R{j} * V;
        RPA{j} = V' * frecov{j} * V;

        %create FRE vec.
        frevec(1,1) = RPA{j}(1,1);
        frevec(2,1) = RPA{j}(1,2);
        frevec(3,1) = RPA{j}(1,3);
        frevec(4,1) = RPA{j}(2,2);
        frevec(5,1) = RPA{j}(2,3);
        frevec(6,1) = RPA{j}(3,3);

        %try moving average.
        %update the mean.
        frecovmean{j} = frecovmean{j} + 1/winsize* (frevec' - frecovwin{j}(end,:));
        %update the window.
        frecovwin{j} = circshift(frecovwin{j},1);
        frecovwin{j}(1,:) = frevec';
        % store for plotting
        frenoisycov{j}(i,:) = frevec;
        fresmoothcov{j}(i,:) = frecovmean{j}';
        fresmoothRMS(i,j) = sqrt(frecovmean{j}(1) + frecovmean{j}(4) + frecovmean{j}(6));
        %estimate the FLE
        %flevec = Ainv{j}*frevec;
        flevec = Ainv{j}*frecovmean{j}';

        %build FLE covariance matrix
        fleest(1,1) = flevec(1,1);
        fleest(1,2) = flevec(2,1);
        fleest(2,1) = flevec(2,1);
        fleest(1,3) = flevec(3,1);
        fleest(3,1) = flevec(3,1);
        fleest(2,2) = flevec(4,1);
        fleest(2,3) = flevec(5,1);
        fleest(3,2) = flevec(5,1);
        fleest(3,3) = flevec(6,1);
        meanFLE = meanFLE + fleest;

        %TODO:rotate back into global axes.
        sigmaFLE = V * fleest * V';

        measfreRMS(i,j) = sqrt(fre(j,:)*fre(j,:)');
        %freRMS(i) = sqrt(curfre'*curfre);
        fleValues{j}(1,i) = fle(1,1);
        fleValues{j}(2,i) = fle(1,2);
        fleValues{j}(3,i) = fle(1,3);
        fleValues{j}(4,i) = fle(2,2);
        fleValues{j}(5,i) = fle(2,3);
        fleValues{j}(6,i) = fle(3,3);

        fleEstValues{j}(1,i) = sigmaFLE(1,1);
        fleEstValues{j}(2,i) = sigmaFLE(1,2);
        fleEstValues{j}(3,i) = sigmaFLE(1,3);
        fleEstValues{j}(4,i) = sigmaFLE(2,2);
        fleEstValues{j}(5,i) = sigmaFLE(2,3);
        fleEstValues{j}(6,i) = sigmaFLE(3,3);

        if(flemodel.heteroscedastic)
            fleEigValues{j}(:,i) = sort(eig(fle(:,:,j)));
            fleRMS{j}(i) = sqrt(trace(fle(:,:,j)));
            meanFLEEigValues(:,i) = sort(eig(fle(:,:,j)));
        else
            fleEigValues{j}(:,i) = sort(eig(fle));
            fleRMS{j}(i) = sqrt(trace(fle));
            meanFLEEigValues(:,i) = sort(eig(fle));
        end
        fleestRMS(i,j) = sqrt(trace(sigmaFLE));
        fleEstEigValues{j}(:,i) = sort(eig(sigmaFLE));
        frecovRMS(i,j) = sqrt(trace(frecov{j}));

    end
    meanFLE = meanFLE/nMrks;
    meanFLERMS(i) = sqrt(trace(meanFLE));
end

fle
timeperframe = (cputime-starttime)/N;
fprintf(fid, 'Time per Frame: %f ms or %f Hz\n', timeperframe*1000, 1/timeperframe);

%% compute the stats.
stats = cell(nMrks);
data = cell(nMrks);
for i = 1:nMrks
    % note we ignore the first set of frames where the covariance matrix is
    % being updated.

    stats{i}.winsize = winsize;
    %RMS.
    stats{i}.RMS.true = fleRMS{i}(1:end);
    stats{i}.RMS.estimated = fleestRMS(1:end,i);
    stats{i}.RMS.error = fleestRMS((winsize+1):end,i) - fleRMS{i}((winsize+1):end);
    stats{i}.RMS.mean = mean(stats{i}.RMS.error);
    stats{i}.RMS.std = std(stats{i}.RMS.error);
    stats{i}.RMS.rms = sqrt(mean((stats{i}.RMS.error).^2));
    stats{i}.RMS.perc95 = [getPercentile(stats{i}.RMS.error,0.025), ...
        getPercentile(stats{i}.RMS.error,0.975)];
    fprintf(fid,'Stats Mrk #%2d\n\tMean\t\tStDev\t\tRMS\t\t95%%Low\t\t95%%Up\nRMS\t% f\t% f\t% f\t% f\t% f\n', i,...
        stats{i}.RMS.mean, stats{i}.RMS.std, stats{i}.RMS.rms, ...
        stats{i}.RMS.perc95(1), stats{i}.RMS.perc95(2));
    fprintf(fid2,'Stats Mrk #%2d\n\tMean\t\tRMS\t\t95%%Low\t\t95%%Up\nRMS\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n', i,...
        stats{i}.RMS.mean, stats{i}.RMS.rms, ...
        stats{i}.RMS.perc95(1), stats{i}.RMS.perc95(2));
    %sigma_11
    stats{i}.sigma11.true = fleValues{i}(1,:);
    stats{i}.sigma11.estimated = fleEstValues{i}(1,:);
    stats{i}.sigma11.error = fleEstValues{i}(1,(winsize+1):end)...
        - fleValues{i}(1,(winsize+1):end);
    stats{i}.sigma11.mean = mean(stats{i}.sigma11.error);
    stats{i}.sigma11.std = std(stats{i}.sigma11.error);
    stats{i}.sigma11.rms = sqrt(mean((stats{i}.sigma11.error).^2));
    stats{i}.sigma11.perc95 = [getPercentile(stats{i}.sigma11.error,0.025),...
        getPercentile(stats{i}.sigma11.error,0.975)];
    fprintf(fid,'S11\t% f\t% f\t% f\t% f\t% f\n',...
        stats{i}.sigma11.mean, stats{i}.sigma11.std, stats{i}.sigma11.rms, ...
        stats{i}.sigma11.perc95(1), stats{i}.sigma11.perc95(2));
    fprintf(fid2,'S11\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n',...
        stats{i}.sigma11.mean, stats{i}.sigma11.rms, ...
        stats{i}.sigma11.perc95(1), stats{i}.sigma11.perc95(2));
    %sigma_12
    stats{i}.sigma12.true = fleValues{i}(2,:);
    stats{i}.sigma12.estimated = fleEstValues{i}(2,:);
    stats{i}.sigma12.error = fleEstValues{i}(2,(winsize+1):end)...
        - fleValues{i}(2,(winsize+1):end);
    stats{i}.sigma12.mean = mean(stats{i}.sigma12.error);
    stats{i}.sigma12.std = std(stats{i}.sigma12.error);
    stats{i}.sigma12.rms = sqrt(mean((stats{i}.sigma12.error).^2));
    stats{i}.sigma12.perc95 = [getPercentile(stats{i}.sigma12.error,0.025),...
        getPercentile(stats{i}.sigma12.error,0.975)];
    fprintf(fid,'S12\t% f\t% f\t% f\t% f\t% f\n',...
        stats{i}.sigma12.mean, stats{i}.sigma12.std, stats{i}.sigma12.rms, ...
        stats{i}.sigma12.perc95(1), stats{i}.sigma12.perc95(2));
    fprintf(fid2,'S12\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n',...
        stats{i}.sigma12.mean, stats{i}.sigma12.rms, ...
        stats{i}.sigma12.perc95(1), stats{i}.sigma12.perc95(2));
    %sigma_13
    stats{i}.sigma13.true = fleValues{i}(3,:);
    stats{i}.sigma13.estimated = fleEstValues{i}(3,:);
    stats{i}.sigma13.error = fleEstValues{i}(3,(winsize+1):end)...
        - fleValues{i}(3,(winsize+1):end);
    stats{i}.sigma13.mean = mean(stats{i}.sigma13.error);
    stats{i}.sigma13.std = std(stats{i}.sigma13.error);
    stats{i}.sigma13.rms = sqrt(mean((stats{i}.sigma13.error).^2));
    stats{i}.sigma13.perc95 = [getPercentile(stats{i}.sigma13.error,0.025),...
        getPercentile(stats{i}.sigma13.error,0.975)];
    fprintf(fid,'S13\t% f\t% f\t% f\t% f\t% f\n',...
        stats{i}.sigma13.mean, stats{i}.sigma13.std, stats{i}.sigma13.rms, ...
        stats{i}.sigma13.perc95(1), stats{i}.sigma13.perc95(2));
    fprintf(fid2,'S13\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n',...
        stats{i}.sigma13.mean, stats{i}.sigma13.rms, ...
        stats{i}.sigma13.perc95(1), stats{i}.sigma13.perc95(2));
    %sigma_22
    stats{i}.sigma22.true = fleValues{i}(4,:);
    stats{i}.sigma22.estimated = fleEstValues{i}(4,:);
    stats{i}.sigma22.error = fleEstValues{i}(4,(winsize+1):end)...
        - fleValues{i}(4,(winsize+1):end);
    stats{i}.sigma22.mean = mean(stats{i}.sigma22.error);
    stats{i}.sigma22.std = std(stats{i}.sigma22.error);
    stats{i}.sigma22.rms = sqrt(mean((stats{i}.sigma22.error).^2));
    stats{i}.sigma22.perc95 = [getPercentile(stats{i}.sigma22.error,0.025),...
        getPercentile(stats{i}.sigma22.error,0.975)];
    fprintf(fid,'S22\t% f\t% f\t% f\t% f\t% f\n',...
        stats{i}.sigma22.mean, stats{i}.sigma22.std, stats{i}.sigma22.rms, ...
        stats{i}.sigma22.perc95(1), stats{i}.sigma22.perc95(2));
    fprintf(fid2,'S22\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n',...
        stats{i}.sigma22.mean, stats{i}.sigma22.rms, ...
        stats{i}.sigma22.perc95(1), stats{i}.sigma22.perc95(2));
    %sigma_23
    stats{i}.sigma23.true = fleValues{i}(5,:);
    stats{i}.sigma23.estimated = fleEstValues{i}(5,:);
    stats{i}.sigma23.error = fleEstValues{i}(5,(winsize+1):end)...
        - fleValues{i}(5,(winsize+1):end);
    stats{i}.sigma23.mean = mean(stats{i}.sigma23.error);
    stats{i}.sigma23.std = std(stats{i}.sigma23.error);
    stats{i}.sigma23.rms = sqrt(mean((stats{i}.sigma23.error).^2));
    stats{i}.sigma23.perc95 = [getPercentile(stats{i}.sigma23.error,0.025),...
        getPercentile(stats{i}.sigma23.error,0.975)];
    fprintf(fid,'S23\t% f\t% f\t% f\t% f\t% f\n',...
        stats{i}.sigma23.mean, stats{i}.sigma23.std, stats{i}.sigma23.rms, ...
        stats{i}.sigma23.perc95(1), stats{i}.sigma23.perc95(2));
    fprintf(fid2,'S23\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n',...
        stats{i}.sigma23.mean, stats{i}.sigma23.rms, ...
        stats{i}.sigma23.perc95(1), stats{i}.sigma23.perc95(2));
    %sigma_33
    stats{i}.sigma33.true = fleValues{i}(6,:);
    stats{i}.sigma33.estimated = fleEstValues{i}(6,:);
    stats{i}.sigma33.error = fleEstValues{i}(6,(winsize+1):end)...
        - fleValues{i}(6,(winsize+1):end);
    stats{i}.sigma33.mean = mean(stats{i}.sigma33.error);
    stats{i}.sigma33.std = std(stats{i}.sigma33.error);
    stats{i}.sigma33.rms = sqrt(mean((stats{i}.sigma33.error).^2));
    stats{i}.sigma33.perc95 = [getPercentile(stats{i}.sigma33.error,0.025),...
        getPercentile(stats{i}.sigma33.error,0.975)];
    fprintf(fid,'S33\t% f\t% f\t% f\t% f\t% f\n',...
        stats{i}.sigma33.mean, stats{i}.sigma33.std, stats{i}.sigma33.rms, ...
        stats{i}.sigma33.perc95(1), stats{i}.sigma33.perc95(2));
    fprintf(fid2,'S33\t% 4.4f & % 4.4f & % 4.4f & % 4.4f\n',...
        stats{i}.sigma33.mean, stats{i}.sigma33.rms, ...
        stats{i}.sigma33.perc95(1), stats{i}.sigma33.perc95(2));
    % stack up the data to be returned for over all analysis.
    data{i}.true = [fleRMS{i}(1:end)';fleValues{i}(1,:);fleValues{i}(2,:);...
        fleValues{i}(3,:); fleValues{i}(4,:); fleValues{i}(5,:); fleValues{i}(6,:)];
    data{i}.estimated = [fleestRMS(1:end,i)';fleEstValues{i}(1,:);fleEstValues{i}(2,:);...
        fleEstValues{i}(3,:); fleEstValues{i}(4,:); fleEstValues{i}(5,:); fleEstValues{i}(6,:)];
end

save(casename, 'stats');

fclose('all');
%% plot the results.
%frenoisycov{1}(6000,:)
if(~bPlot)
    return;
end

figure(1);
clf;
if(flemodel.heteroscedastic)
    plot(1:N, fleRMS{1}, '-', 'Color', (2/3)*[1,1,1], 'LineWidth',2);
    hold on;
    plot(1:N, fleRMS{2}, '-', 'Color', (2/3)*[1,1,1], 'LineWidth',1);
    plot(1:N, fleRMS{3}, '-.', 'Color', (2/3)*[1,1,1], 'LineWidth',2);
    plot(1:N, fleRMS{4}, '--', 'Color', (2/3)*[1,1,1], 'LineWidth',2);
    if(nMrks > 4)
        plot(1:N, fleRMS{5}, 'k-', 'LineWidth',2);
    end
else
    plot(1:N, fleRMS{1}, 'k:', 'LineWidth',2);
    hold on;
end
% plot(1:N, frecovRMS(:,1), 'r');
% plot(1:N, frecovRMS(:,2), 'g');
% plot(1:N, frecovRMS(:,3), 'b');
% plot(1:N, frecovRMS(:,4), 'y');
% if(nMrks > 4)
%     plot(1:N, frecovRMS(:,5), 'c');
% end
hold off;
titlestring = sprintf('Heteroscedastic FLE');
title(titlestring, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('RMS (mm)', 'fontsize', 18);
if(nMrks > 4)
    legend('Mrk #1', 'Mrk #2', 'Mrk #3', 'Mrk #4', 'Mrk #5', 'Location', 'SouthEast');
else
    legend('Mrk #1', 'Mrk #2', 'Mrk #3', 'Mrk #4', 'Location', 'SouthEast');
end
set(gca, 'fontsize', 14);
filename = sprintf('%s-trueFLE-AllMarkers.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-trueFLE-AllMarkers.jpg', casename);
print('-djpeg', filename);


figure(2);
clf;
if(flemodel.heteroscedastic)
    plot(1:N, fleRMS{1}, '-', 'Color', (2/3)*[1,1,1], 'LineWidth',2);
    hold on;
    plot(1:N, fleRMS{2}, '-', 'Color', (2/3)*[1,1,1], 'LineWidth',1);
    plot(1:N, fleRMS{3}, '-.', 'Color', (2/3)*[1,1,1], 'LineWidth',2);
    plot(1:N, fleRMS{4}, '--', 'Color', (2/3)*[1,1,1], 'LineWidth',2);
    if(nMrks > 4)
        plot(1:N, fleRMS{5}, 'c:', 'LineWidth', 2);
    end
else
    plot(1:N, fleRMS{1}, 'k--', 'LineWidth',2);
    hold on;
end
plot(1:N, fleestRMS(:,1), '-', 'Color', (1/3)*[1, 1, 1], 'LineWidth',2);
plot(1:N, fleestRMS(:,2), '-', 'Color', (1/3)*[1, 1, 1], 'LineWidth',1);
plot(1:N, fleestRMS(:,3), '-.', 'Color', (1/3)*[1, 1, 1], 'LineWidth',2);
plot(1:N, fleestRMS(:,4), '--', 'Color', (1/3)*[1, 1, 1], 'LineWidth',2);
if(nMrks > 4)
    plot(1:N, fleestRMS(:,5), '-', 'Color', (2/3)*[1, 1, 1], 'LineWidth',2);
end
plot(1:N, meanFLERMS, 'k', 'LineWidth', 2);
%plot(1:N, measfreRMS, '.b');
%plot(1:N, sqrt(2)*RRMS(:,1), 'c');
hold off;
ylim([0 0.55]);
titlestring1 = sprintf('Estimated FLE RMS');
titlestring = sprintf('(FRE computed with %d frames window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('RMS (mm)', 'fontsize', 18);
if(flemodel.heteroscedastic)
    if(nMrks > 4)
        legend('True FLE Mrk #1', 'True FLE Mrk #2', 'True FLE Mrk #3',...
            'True FLE Mrk #4', 'True FLE Mrk #5', 'Est. FLE Mrk #1',...
            'Est. FLE Mrk #2', 'Est. FLE Mrk #3', 'Est. FLE Mrk #4',...
            'Est. FLE Mrk #5', 'Mean of Markers', 'Location', 'SouthEast');
    else
        legend('True FLE Mrk #1', 'True FLE Mrk #2', 'True FLE Mrk #3',...
            'True FLE Mrk #4', 'Est. FLE Mrk #1', 'Est. FLE Mrk #2',...
            'Est. FLE Mrk #3', 'Est. FLE Mrk #4',...
            'Mean of Markers', 'Location', 'SouthEast');
    end
else
    if(nMrks > 4)
        legend('True FLE RMS', 'Mrk #1', 'Mrk #2', 'Mrk #3', 'Mrk #4', 'Mrk #5',...
            'Mean of Markers', 'Location', 'SouthEast');
    else
        legend('True FLE RMS', 'Mrk #1', 'Mrk #2', 'Mrk #3', 'Mrk #4', ...
            'Mean of Markers', 'Location', 'SouthEast');
    end
end
set(gca, 'fontsize', 14);
filename = sprintf('%s-estimatedFLE-AllMarkers.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLE-AllMarkers.jpg', casename);
print('-djpeg',filename);

figure(3); % show the computation on Marker 1.
plot(1:N, fleRMS{1}, 'k--', 'LineWidth', 2);
hold on;
plot(1:N, frecovRMS(:,1), '-', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 1);
plot(1:N, fresmoothRMS(:,1), '-.', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 2);
plot(1:N, fleestRMS(:,1), 'k', 'LineWidth', 2);
%plot(1:N, meanFLERMS, 'c');
%plot(1:N, measfreRMS, '.b');
plot(1:N, sqrt(nMrks/(nMrks-2))*fresmoothRMS(:,1), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.55]);
titlestring1 = sprintf('Windowed FRE and Estimated FLE for Marker #1');
titlestring = sprintf('(FRE computed with %d frames window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('RMS (mm)', 'fontsize', 18);
legend('True FLE RMS', 'Windowed FRE', 'Moving Average of FRE', 'Estimated FLE',...
    'Estimated FLE using Fitzpatrick 1998', 'Location', 'SouthEast');
set(gca, 'fontsize', 14);
filename = sprintf('%s-estimatedFLE-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLE-Mrk01.jpg', casename);
print('-djpeg',filename);

figure(4); % plot the eigenvalues.
plot(1:N, fleEigValues{1}(1,:), 'r:');
hold on;
plot(1:N, fleEstEigValues{1}(1,:), 'r');
plot(1:N, fleEigValues{1}(2,:), 'g:');
plot(1:N, fleEstEigValues{1}(2,:), 'g');
plot(1:N, fleEigValues{1}(3,:), 'b:');
plot(1:N, fleEstEigValues{1}(3,:), 'b');
hold off;

figure(5);
plot(1:N, frecovRMS(:,1), 'b');
hold on;
plot(1:N, fresmoothRMS(:,1), 'g');
hold off;

figure(6);
for i = 1:6
    subplot(6,1,i);
    plot(1:N, frenoisycov{1}(:,i), 'b');
    hold on;
    plot(1:N, fresmoothcov{1}(:,i), 'g', 'LineWidth', 2);
    hold off;
end

figure(7);
plot(1:N, fleValues{1}(1,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(1,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.25]);
title('Estimated FLE Component: \sigma_{11}', 'fontsize', 20);
set(gca, 'fontsize', 18);
filename = sprintf('%s-estimatedFLESigma11-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLESigma11-Mrk01.jpg', casename);
print('-djpeg',filename);

figure(8);
plot(1:N, fleValues{1}(2,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(2,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([-0.025 0.025]);
title('Estimated FLE Component: \sigma_{12}', 'fontsize', 20);
set(gca, 'fontsize', 18);
filename = sprintf('%s-estimatedFLESigma12-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLESigma12-Mrk01.jpg', casename);
print('-djpeg',filename);

figure(9);
plot(1:N, fleValues{1}(3,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(3,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([-0.025 0.025]);
title('Estimated FLE Component: \sigma_{13}', 'fontsize', 20);
set(gca, 'fontsize', 18);
filename = sprintf('%s-estimatedFLESigma13-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLESigma13-Mrk01.jpg', casename);
print('-djpeg',filename);

figure(10);
plot(1:N, fleValues{1}(4,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(4,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.25]);
title('Estimated FLE Component: \sigma_{22}', 'fontsize', 20);
set(gca, 'fontsize', 18);
filename = sprintf('%s-estimatedFLESigma22-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLESigma22-Mrk01.jpg', casename);
print('-djpeg',filename);

figure(11);
plot(1:N, fleValues{1}(5,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(5,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([-0.025 0.025]);
title('Estimated FLE Component: \sigma_{23}', 'fontsize', 20);
set(gca, 'fontsize', 18);
filename = sprintf('%s-estimatedFLESigma23-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLESigma23-Mrk01.jpg', casename);
print('-djpeg',filename);

figure(12);
plot(1:N, fleValues{1}(6,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(6,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.25]);
title('Estimated FLE Component: \sigma_{33}', 'fontsize', 20);
set(gca, 'fontsize', 18);
filename = sprintf('%s-estimatedFLESigma33-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-estimatedFLESigma33-Mrk01.jpg', casename);
print('-djpeg',filename);
