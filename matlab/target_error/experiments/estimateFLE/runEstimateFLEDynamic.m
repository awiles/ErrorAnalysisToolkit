%runEstimateFLEDynamic
email_setup;
%cd E:\docs\research\phd\experiments\FLEPrediction\simulated-tool-paths
cd E:\awiles\data\FLEPrediction
%% 1. Tool design.
%body = 'west';
%body = 'west2';
body = 'ta003-4';
%get the rigid body design.
switch(body)
    case 'west'
        casename = 'WestTool';
        [refmrk, normals, tip] = getWestToolDesign('e', 71, 54, 85);
    case 'west2'
        casename = 'WestTool2';
        [refmrk, normals, tip] = getWestToolDesign('e', 71, 54, 85);
        refmrk = [refmrk; 0, 0, 5];
    case 'ta003-4'
        casename = 'TA003-4';
        refmrk = [0.0000, 0.0000, -0.0120;...
            -99.4066, -33.4878, 0.0417;...
            -153.6516, 0.2701, -0.0463;...
            -103.4595, 49.1463, 0.0287];
        %         refmrk = [0.0000, 0.0000, 0.0;...
        %             -99.4066, -33.4878, 0.0;...
        %             -153.6516, 0.2701,  0.0;...
        %             -103.4595, 49.1463, 0.0];
        tip = [200, 0, 0];
    otherwise
        error('Invalid body design.');
end
tooldesign.name = casename;
tooldesign.refmrk = refmrk;
tooldesign.tip = tip;
[U, L, tooldesign.V0, x, tooldesign.A, tooldesign.Ainv] = xfrmToPA(refmrk);

%% 2. FLE Model.
model = 'quad1';
switch(model)
    case 'quad1'
        flemodel.name = 'RadialQuadratic500um';
        flemodel.model = 'radial';
        flemodel.rate = [5.7e-9, 5.7e-9, 5.1e-8];
    otherwise
        error('Invalide FLE model given.');
end
%% 3. Tool Path.
path = 3;
filename = sprintf('path-%02d', path);
load(filename);
toolpath.name = filename;
toolpath.xfrm = xfrm;

%% 4. Window Size.
winsize = 500;

%% 5. Number of samples.
M = 500;

%% set up the simulation.
N = length(xfrm);
sumdata = zeros(7, N, M);
casename = sprintf('%s-%s-%s-%03d', tooldesign.name, flemodel.name,...
    toolpath.name, winsize);
mkdir(casename);
cd(casename);

%% Save the parameters.
parmfilename = sprintf('%s-parameters', casename);
save(parmfilename, 'tooldesign', 'flemodel', 'toolpath', 'winsize');

%% run the simulation M times.
starttime = cputime;
for i = 1:M
    data = estimateFLEDynamic(tooldesign, flemodel, toolpath, winsize, i, 1);
    sumdata(:, :, i) = [data{1}.estimated];
    if(i==1)
        truedata = data{1}.true;
    end
    
    % compute estimated completion time.
    avgTime = (cputime - starttime)/i; %in seconds.
    estTimeToFinish = avgTime*(M-i); %in seconds.
    estHours = floor(estTimeToFinish/3600);
    estMinutes = floor((estTimeToFinish - 3600*estHours)/60);
   
    fprintf('Completed %3.1f %%... Avg. time per run: %4.2f seconds. Finish in approx. %d hours %d min \n',...
        (100*i/M), avgTime, estHours, estMinutes);
    if(rem(i,100) == 0 || i == 1)
        msg = sprintf('Completed %3.1f percent... Avg. time per run: %4.2f seconds. Finish in approx. %d hours %d min \n',...
        (100*i/M), avgTime, estHours, estMinutes);
        curlogfile = sprintf('%s-%04d.txt', casename, i);
        curRMSplot = sprintf('%s-%04d-estimatedFLE-AllMarkers.eps', casename, i);
        curRMSplot2 = sprintf('%s-%04d-estimatedFLE-Mrk01.eps', casename, i);
        sendmail('adwiles@gmail.com', 'FLE Simulation Update', msg, {curlogfile, curRMSplot, curRMSplot2});
    end

end

save(casename, 'sumdata', 'truedata');
clear sumdata truedata;
[logfilename, filenamerms] = plotEstimateFLEDynamicSummaryData(casename, {'Best', 'Best'});
sendmail('adwiles@gmail.com', 'FLE Simulation Complete', 'FLE Simulation Complete', {logfilename,filenamerms});

cd ..

