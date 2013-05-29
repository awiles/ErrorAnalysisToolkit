function plotRandFiducialTargetExperiment(testname, datetime, fleModel, figFontSize, varargin)

nToPlot = -1;
verbose = 0;

if( nargin > 4 )
    nVarArgs = length(varargin);
    i = 1;
    while( i <= nVarArgs )
        if( strcmp(varargin{i}, 'PlotDataSamples') )
            i=i+1;
            nToPlot = varargin{i};
            if( ~isscalar(nToPlot) && ~isvector(nToPlot) )
                error('Value for PlotDataSamples is invalid');
            end
        elseif( strcmp(varargin{i}, 'Verbose'))
            verbose = 1;
        else
            error('Unknown paramter: %s', varargin{i});
        end
        i=i+1;
    end
end

cd(testname)
cd(datetime)

data = csvread('data.csv');
id = 1;
covonly = 2;
muandcov = 3;
measRMS = 4;
theoryRMS = 5;
RMSDiff = 6;
percentDiff = 7;

load parm;

plot(data(:,theoryRMS), data(:,measRMS), '.k');
hold on;
minRMS = min(min(data(:,measRMS:theoryRMS)));
maxRMS = max(max(data(:,measRMS:theoryRMS)));
plot([minRMS maxRMS], [minRMS,maxRMS], 'k-');
hold off;

xlabel('Predicted TRE RMS (mm)', 'fontsize',figFontSize);
ylabel('Simulated TRE RMS (mm)', 'fontsize',figFontSize);
set(gca,'fontsize',figFontSize);

titlestring = sprintf('Predicted versus Simulated TRE RMS of Random Configurations');
subtitlestring = sprintf('%s FLE Model, RMS_{FLE} = %3.2f, %d Fiducials',...
    fleModel, sqrt(trace(parm.Sigma)), parm.nMarkers);
title({titlestring, subtitlestring}, 'fontsize',figFontSize);

axis([minRMS, maxRMS, minRMS, maxRMS]);
figurefilename = sprintf('PredvMeas_%s', parm.name);
%print('-depsc', '-tiff', '-r300', figurefilename);
print('-dpng', '-r600', figurefilename);

%copyfile('*.eps', 'E:\awiles\data\tretest\IEEE_Data\RandDesigns' );

%% lets plot the histograms for the distance error for each case.
nTotalCount = parm.nBodies*parm.nTrials*parm.nOrientations*parm.nPositions;
%histbins = 0.01:0.01:1

if( isscalar(nToPlot) && nToPlot < 1 )
    nToPlot = 1:nTotalCount;
end

for i = nToPlot
    filename = sprintf('data%06d.mat', i);
    fprintf('Plotting for %s ...\n', filename);
    load(filename);
    % generate histogram for theoterical data.
    
    % generate experimental histogram.
    disterr = sqrt(sum((probe.Meas.error).^2,2));
    %mchistdata = hist(disterr, histbins);
    % plot each one.
    figure(2);    
    [mu, sigma, rms, histbins] = plotCumHist(disterr,[], '.-k');
    hold on;
    %plotMaxwellDist( testResult.probe.RMS.Meas, histbins, '-r');
    plotRadius3DDist( probe.Meas.stats.Sigma, histbins, '-r');
    probe.Meas.stats.Sigma;
    %plotMaxwellDist( testResult.probe.RMS.Theory, histbins, '-g');
    plotRadius3DDist( probe.Actual.stats.Sigma, histbins, '-g');
    probe.Actual.stats.Sigma;
    %plotMaxwellDist( testResult.probeAvgFLE.RMS.Theory, histbins, '-b');
    %plotRadius3DDist( probe.AvgFLE.stats.Sigma, histbins, '-b');
    %probe.AvgFLE.stats.Sigma;
    hold off;
    xlim([0, max(histbins)]);
    legend('Monte Carlo', 'Measured TRE', 'NH TRE');
    titlestring = sprintf('TRE Distribution %06d', i);
    subtitlestring = sprintf('Meas. RMS: %3.2f, Theory RMS: %3.2f', ...
        testResult.probe.RMS.Meas, testResult.probe.RMS.Theory);
    title({titlestring, subtitlestring});
    figurefilename = sprintf('hist_%s_%06d', parm.name, i);
    %print('-depsc', '-tiff', '-r300', figurefilename);
    print('-dpng', '-r600', figurefilename);
end

cd ../..
