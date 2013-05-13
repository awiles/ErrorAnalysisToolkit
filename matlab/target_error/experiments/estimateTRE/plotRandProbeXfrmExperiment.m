function plotRandFiducialTargetExperiment(testname, fleModel, figFontSize)

cd(testname)

data = csvread('data.csv');
id = 1;
covonly = 2;
muandcov = 3;
measRMS = 4;
theoryRMS = 5;
RMSDiff = 6;
percentDiff = 7;



load parm;

%% histogram
figure(1);
hist(data(:,percentDiff));
axis([-1,1,0,300]);
titlestring = sprintf('Histogram of the %% Differences for RMS');
subtitlestring = sprintf('%s FLE Model, RMS_{FLE} = %3.2f, Tool %s, {\\rho} = %d mm',...
    fleModel, sqrt(trace(parm.Sigma)), parm.design, parm.rho);
title({titlestring, subtitlestring}, 'fontsize',figFontSize);
set(gca,'fontsize',figFontSize);
figurefilename = sprintf('Histogram_%s', parm.name);
print('-depsc', '-tiff', '-r300', figurefilename);

%% plot the scatter plot.
figure(2);
plot(data(:,theoryRMS), data(:,measRMS), '.k');
hold on;
minRMS = min(min(data(:,measRMS:theoryRMS)));
maxRMS = max(max(data(:,measRMS:theoryRMS)));
plot([minRMS maxRMS], [minRMS,maxRMS], 'k-');
hold off;

xlabel('Predicted TRE RMS (mm)', 'fontsize',figFontSize);
ylabel('Simulated TRE RMS (mm)', 'fontsize',figFontSize);
set(gca,'fontsize',figFontSize);

titlestring = sprintf('Predicted versus Simulated TRE RMS');
subtitlestring = sprintf('%s FLE Model, RMS_{FLE} = %3.2f, Tool %s, {\\rho} = %d mm',...
    fleModel, sqrt(trace(parm.Sigma)), parm.design, parm.rho);
title({titlestring, subtitlestring}, 'fontsize',figFontSize);

axis([minRMS, maxRMS, minRMS, maxRMS]);
figurefilename = sprintf('PredvMeas_%s', parm.name);
print('-depsc', '-tiff', '-r300', figurefilename);

copyfile('*.eps', 'E:\awiles\data\tretest\IEEE_Data\WestTools' );

cd ..
