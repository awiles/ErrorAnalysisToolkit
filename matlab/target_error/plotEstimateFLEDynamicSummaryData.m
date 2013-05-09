function [logfilename, filenamerms] = plotEstimateFLEDynamicSummaryData(casename, legPos, bounds)
%plot the summary data.
load(casename);
parmfilename = sprintf('%s-parameters', casename);
load(parmfilename);

N = size(sumdata,2);
M = size(sumdata,3);

%% post process the summary data
rmsdata = zeros(3,N);
sigma11data = zeros(3,N);
sigma12data = zeros(3,N);
sigma13data = zeros(3,N);
sigma22data = zeros(3,N);
sigma23data = zeros(3,N);
sigma33data = zeros(3,N);
for i = 1:N
    %rms
    rmsmean = mean(sumdata(1,i,:));
    rmsperc025 = getPercentile(sumdata(1,i,:),0.025);
    rmsperc975 = getPercentile(sumdata(1,i,:),0.975);
    rmsdata(:,i) = [rmsmean; rmsperc025; rmsperc975];
    %sigma11
    sigma11mean = mean(sumdata(2,i,:));
    sigma11perc025 = getPercentile(sumdata(2,i,:),0.025);
    sigma11perc975 = getPercentile(sumdata(2,i,:),0.975);
    sigma11data(:,i) = [sigma11mean; sigma11perc025; sigma11perc975];
    %sigma12
    sigma12mean = mean(sumdata(3,i,:));
    sigma12perc025 = getPercentile(sumdata(3,i,:),0.025);
    sigma12perc975 = getPercentile(sumdata(3,i,:),0.975);
    sigma12data(:,i) = [sigma12mean; sigma12perc025; sigma12perc975];
    %sigma13
    sigma13mean = mean(sumdata(4,i,:));
    sigma13perc025 = getPercentile(sumdata(4,i,:),0.025);
    sigma13perc975 = getPercentile(sumdata(4,i,:),0.975);
    sigma13data(:,i) = [sigma13mean; sigma13perc025; sigma13perc975];
    %sigma22
    sigma22mean = mean(sumdata(5,i,:));
    sigma22perc025 = getPercentile(sumdata(5,i,:),0.025);
    sigma22perc975 = getPercentile(sumdata(5,i,:),0.975);
    sigma22data(:,i) = [sigma22mean; sigma22perc025; sigma22perc975];
    %sigma23
    sigma23mean = mean(sumdata(6,i,:));
    sigma23perc025 = getPercentile(sumdata(6,i,:),0.025);
    sigma23perc975 = getPercentile(sumdata(6,i,:),0.975);
    sigma23data(:,i) = [sigma23mean; sigma23perc025; sigma23perc975];
    %sigma33
    sigma33mean = mean(sumdata(7,i,:));
    sigma33perc025 = getPercentile(sumdata(7,i,:),0.025);
    sigma33perc975 = getPercentile(sumdata(7,i,:),0.975);
    sigma33data(:,i) = [sigma33mean; sigma33perc025; sigma33perc975];
end

%% compute span stats.
logfilename = sprintf('%s.txt',casename);
fid = fopen(logfilename, 'wt');
fprintf(fid,'Average Statistics Range\n');
fprintf(fid,'\tSpan: Upper 95%% Bound less the Lower 95%% Bound\n');
fprintf(fid,'\tAvg. Up: Upper 95%% Bound less Mean\n');
fprintf(fid,'\tAvg. Low: Mean less the Lower 95%% Bound\n\n');
fprintf(fid,'\tSpan\t\tAvg. Up\t\tAvg. Low\n');
avgSpanRMS = mean(rmsdata(3,(2*winsize+1):end) - rmsdata(2,(2*winsize+1):end));
avgUpperRMS = mean(rmsdata(3,(2*winsize+1):end) - rmsdata(1,(2*winsize+1):end));
avgLowerRMS = mean(rmsdata(2,(2*winsize+1):end) - rmsdata(1,(2*winsize+1):end));
fprintf(fid,'RMS:\t% f\t% f\t% f\n', avgSpanRMS, avgUpperRMS, avgLowerRMS);
avgSpanSigma11 = mean(sigma11data(3,(2*winsize+1):end) - sigma11data(2,(2*winsize+1):end));
avgUpperSigma11 = mean(sigma11data(3,(2*winsize+1):end) - sigma11data(1,(2*winsize+1):end));
avgLowerSigma11 = mean(sigma11data(2,(2*winsize+1):end) - sigma11data(1,(2*winsize+1):end));
fprintf(fid,'S11:\t% f\t% f\t% f\n', avgSpanSigma11, avgUpperSigma11, avgLowerSigma11);
avgSpanSigma12 = mean(sigma12data(3,(2*winsize+1):end) - sigma12data(2,(2*winsize+1):end));
avgUpperSigma12 = mean(sigma12data(3,(2*winsize+1):end) - sigma12data(1,(2*winsize+1):end));
avgLowerSigma12 = mean(sigma12data(2,(2*winsize+1):end) - sigma12data(1,(2*winsize+1):end));
fprintf(fid,'S12:\t% f\t% f\t% f\n', avgSpanSigma12, avgUpperSigma12, avgLowerSigma12);
avgSpanSigma13 = mean(sigma13data(3,(2*winsize+1):end) - sigma13data(2,(2*winsize+1):end));
avgUpperSigma13 = mean(sigma13data(3,(2*winsize+1):end) - sigma13data(1,(2*winsize+1):end));
avgLowerSigma13 = mean(sigma13data(2,(2*winsize+1):end) - sigma13data(1,(2*winsize+1):end));
fprintf(fid,'S13:\t% f\t% f\t% f\n', avgSpanSigma13, avgUpperSigma13, avgLowerSigma13);
avgSpanSigma22 = mean(sigma22data(3,(2*winsize+1):end) - sigma22data(2,(2*winsize+1):end));
avgUpperSigma22 = mean(sigma22data(3,(2*winsize+1):end) - sigma22data(1,(2*winsize+1):end));
avgLowerSigma22 = mean(sigma22data(2,(2*winsize+1):end) - sigma22data(1,(2*winsize+1):end));
fprintf(fid,'S22:\t% f\t% f\t% f\n', avgSpanSigma22, avgUpperSigma22, avgLowerSigma22);
avgSpanSigma23 = mean(sigma23data(3,(2*winsize+1):end) - sigma23data(2,(2*winsize+1):end));
avgUpperSigma23 = mean(sigma23data(3,(2*winsize+1):end) - sigma23data(1,(2*winsize+1):end));
avgLowerSigma23 = mean(sigma23data(2,(2*winsize+1):end) - sigma23data(1,(2*winsize+1):end));
fprintf(fid,'S23:\t% f\t% f\t% f\n', avgSpanSigma23, avgUpperSigma23, avgLowerSigma23);
avgSpanSigma33 = mean(sigma33data(3,(2*winsize+1):end) - sigma33data(2,(2*winsize+1):end));
avgUpperSigma33 = mean(sigma33data(3,(2*winsize+1):end) - sigma33data(1,(2*winsize+1):end));
avgLowerSigma33 = mean(sigma33data(2,(2*winsize+1):end) - sigma33data(1,(2*winsize+1):end));
fprintf(fid,'S33:\t% f\t% f\t% f\n', avgSpanSigma33, avgUpperSigma33, avgLowerSigma33);

fprintf(fid, '\nCovariance and RMS at Frame %d\n\tTrue\tMean\t95Low\t95Up\n', N/2);
fprintf(fid, 'RMS:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(1,N/2), rmsdata(1,N/2), rmsdata(2,N/2), rmsdata(3,N/2));
fprintf(fid, 'S11:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(2,N/2), sigma11data(1,N/2), sigma11data(2,N/2), sigma11data(3,N/2));
fprintf(fid, 'S12:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(3,N/2), sigma12data(1,N/2), sigma12data(2,N/2), sigma12data(3,N/2));
fprintf(fid, 'S13:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(4,N/2), sigma13data(1,N/2), sigma13data(2,N/2), sigma13data(3,N/2));
fprintf(fid, 'S22:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(5,N/2), sigma22data(1,N/2), sigma22data(2,N/2), sigma22data(3,N/2));
fprintf(fid, 'S23:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(6,N/2), sigma23data(1,N/2), sigma23data(2,N/2), sigma23data(3,N/2));
fprintf(fid, 'S33:\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n', truedata(7,N/2), sigma33data(1,N/2), sigma33data(2,N/2), sigma33data(3,N/2));


fclose('all');

%% Plot the RMS results.
figure(20);
plot(1:N, truedata(1,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, rmsdata(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, rmsdata(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, rmsdata(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([0 bounds(1)]);
end
titlestring1 = sprintf('Average RMS FLE over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('RMS (mm)', 'fontsize', 18);
legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
    'Location', legPos{1});
set(gca, 'fontsize', 14);
filenamerms = sprintf('%s-averageFLERMS-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filenamerms);
filenamerms = sprintf('%s-averageFLERMS-Mrk01.png', casename);
print('-dpng','-r600', filenamerms);

%% Sigma11.
figure(21);
plot(1:N, truedata(2,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, sigma11data(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, sigma11data(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, sigma11data(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([0 bounds(2)]);
end
titlestring1 = sprintf('Average \\sigma_{11} over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('\sigma_{11} (mm^2)', 'fontsize', 18);
legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
    'Location', legPos{2});
set(gca, 'fontsize', 14);
filename = sprintf('%s-averageFLESigma11-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-averageFLESigma11-Mrk01.png', casename);
print('-dpng','-r600', filename);

%% Sigma12.
figure(22);
plot(1:N, truedata(3,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, sigma12data(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, sigma12data(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, sigma12data(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([-1*bounds(3) bounds(3)]);
end
titlestring1 = sprintf('Average \\sigma_{12} over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('\sigma_{12} (mm^2)', 'fontsize', 18);
%legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
%    'Location', legPos);
set(gca, 'fontsize', 14);
filename = sprintf('%s-averageFLESigma12-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-averageFLESigma12-Mrk01.png', casename);
print('-dpng','-r600', filename);

%% Sigma13.
figure(23);
plot(1:N, truedata(4,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, sigma13data(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, sigma13data(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, sigma13data(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([-1*bounds(3) bounds(3)]);
end
titlestring1 = sprintf('Average \\sigma_{13} over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('\sigma_{13} (mm^2)', 'fontsize', 18);
%legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
%    'Location', legPos);
set(gca, 'fontsize', 14);
filename = sprintf('%s-averageFLESigma13-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-averageFLESigma13-Mrk01.png', casename);
print('-dpng','-r600', filename);


%% Sigma22.
figure(24);
plot(1:N, truedata(5,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, sigma22data(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, sigma22data(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, sigma22data(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([0 bounds(2)]);
end
titlestring1 = sprintf('Average \\sigma_{22} over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('\sigma_{22} (mm^2)', 'fontsize', 18);
%legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
%    'Location', legPos);
set(gca, 'fontsize', 14);
filename = sprintf('%s-averageFLESigma22-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-averageFLESigma22-Mrk01.png', casename);
print('-dpng','-r600', filename);

%% Sigma23.
figure(25);
plot(1:N, truedata(6,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, sigma23data(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, sigma23data(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, sigma23data(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([-1*bounds(3) bounds(3)]);
end
titlestring1 = sprintf('Average \\sigma_{23} over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('\sigma_{23} (mm^2)', 'fontsize', 18);
%legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
%    'Location', legPos);
set(gca, 'fontsize', 14);
filename = sprintf('%s-averageFLESigma23-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-averageFLESigma23-Mrk01.png', casename);
print('-dpng','-r600', filename);

%% Sigma33.
figure(26);
plot(1:N, truedata(7,:), 'k--', 'LineWidth', 2); %true.
hold on;
plot((2*winsize+1):N, sigma33data(1,(2*winsize+1):N), 'k-', 'LineWidth', 2);
plot((2*winsize+1):N, sigma33data(2,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([2/3, 2/3, 2/3]));
plot((2*winsize+1):N, sigma33data(3,(2*winsize+1):N), '-', 'LineWidth', 2, 'Color', ([1/3, 1/3, 1/3]));
hold off;
if(nargin > 2)
    ylim([0 bounds(2)]);
end
titlestring1 = sprintf('Average \\sigma_{33} over %d Simulations', M);
titlestring = sprintf('(Marker #1, %d frames sliding window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('\sigma_{33} (mm^2)', 'fontsize', 18);
%legend('True', 'Average', '95% Lower Bound', '95% Upper Bound',...
%    'Location', legPos);
set(gca, 'fontsize', 14);
filename = sprintf('%s-averageFLESigma33-Mrk01.eps', casename);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('%s-averageFLESigma33-Mrk01.png', casename);
print('-dpng','-r600', filename);