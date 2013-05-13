%runEstimateFLEDynamic
%email_setup;
cd E:\docs\research\phd\experiments\FLEPrediction\patent-plots
%cd E:\awiles\data\FLEPrediction
%% 1. Tool design.
%body = 'west';
body = 'west2';
%body = 'ta003-4';
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
        flemodel.heteroscedastic = 0;
    case 'quad2'
        flemodel.name = 'RadialQuadratic500umHS';
        flemodel.model = 'radial';
        flemodel.rate = [5.7e-9, 5.7e-9, 5.1e-8];
        flemodel.heteroscedastic = 1;
    otherwise
        error('Invalid FLE model given.');
end
%% 3. Tool Path.
path = 1;
filename = sprintf('path-%02d', path);
load(filename);
toolpath.name = filename;
toolpath.xfrm = xfrm;

%% 4. Window Size.
winsize = 500;

%% set up the simulation.
N = length(xfrm);
casename = sprintf('%s-%s-%s-%03d', tooldesign.name, flemodel.name,...
    toolpath.name, winsize);

%% Save the parameters.
parmfilename = sprintf('%s-parameters', casename);
save(parmfilename, 'tooldesign', 'flemodel', 'toolpath', 'winsize');

%% run the simulation
[data, stats] = estimateFLEDynamic(tooldesign, flemodel, toolpath, winsize, 1001, 1);
stats4 = stats;