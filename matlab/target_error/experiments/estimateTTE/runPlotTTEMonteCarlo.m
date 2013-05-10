function runPlotTTEMonteCarlo(name, test, axisLim)

close all;
cd(name);
currentcasename = sprintf('%s-%04d', name, test);
load(currentcasename);
if(nargin > 2)
    plotTTEMonteCarlo(data,stats,currentcasename, axisLim);
else
    plotTTEMonteCarlo(data,stats,currentcasename);
end
copyfile('*.png', '\\San-350\awiles\data\TTE6DImages');

cd ..