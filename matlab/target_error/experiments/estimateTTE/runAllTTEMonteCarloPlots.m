%% do all plots.
axisLim{1} = [-0.75 0.75 -0.75 0.75 -0.75 0.75];
axisLim{2} = [-0.75 0.75 -0.75 0.75 -0.75 0.75];
axisLim{3} = [-0.75 0.75 -0.75 0.75 -0.75 0.75];
axisLim{4} = [-0.75 0.75 -0.75 0.75 -0.75 0.75];
axisLim{5} = [-3 3 -3 3 -3 3];
cd E:\awiles\data\TTEMonteCarlo
for i=1:3
    for j=1:5
        %casename = sprintf('TTENeedle-tar%02d-xfrm%02d', i, j);
        casename = sprintf('TTERadolfzell-tar%02d-xfrm%02d', i, j);
        runPlotTTEMonteCarlo(casename,1, axisLim{j}); %run the first test for each case.
    end
end