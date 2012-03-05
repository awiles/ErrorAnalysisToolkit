% plot the four covariance matrices (true, mean, 95% upper and lower
% bounds).
close all;
cd E:\docs\research\phd\experiments\FLEPrediction\patent-plots
fsize = 24;
body = 'west';
%body = 'west2';
%body = 'ta003-4';
switch(body)
    case 'west'
        casename = 'WestTool';
        [refmrk, normals, tip] = getWestToolDesign('e', 71, 54, 85);
        S_true = [0.0129, 0.0000, 0.0000; 0.0000, 0.0129, 0.0000; 0.0000, 0.0000, 0.1155];
        S_mean = [0.0130, 0.0001, 0.0002; 0.0001, 0.0130, 0.0017; 0.0002, 0.0017, 0.1163];
        S_95L  = [0.0087, -0.0033, -0.0093; -0.0033, 0.0097, -0.0086; -0.0093, -0.0086, 0.0923];
        S_95U  = [0.0177, 0.0029, 0.0110; 0.0029, 0.0166, 0.0120; 0.0110, 0.0120, 0.1436];
        load path-01
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

%% rotate the markers into the given xfrm

mrk = (xfrm{6000}.R * refmrk')' + repmat(xfrm{6000}.pos, size(refmrk,1), 1);
tip0 = (xfrm{6000}.R * tip')' + repmat(xfrm{6000}.pos, size(tip,1), 1);

conf = 0.95; % confidence regions.

figure(1);
clf;
% YZ -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_true(2:3,2:3), 'k', '--');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_mean(2:3,2:3), 'k', '-');
hold on;
area.C95L = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_95L(2:3,2:3), 2/3*[1,1,1], '-');
hold on;
area.C95U = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_95U(2:3,2:3), 1/3*[1,1,1], '-.');
xlabel('Y (mm)','fontsize', fsize);
ylabel('Z (mm)','fontsize', fsize);
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', fsize);
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeFLEYZPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeFLEYZPlot.jpg');

fprintf('FLE YZ Plane Area Perc. Diff.  Low: %f, Up: %f\n', ...
    100*(area.C95L-area.true)/area.true, 100*(area.C95U-area.true)/area.true);

figure(2);
clf;
% XZ -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [S_true(1,1) S_true(1,3); S_true(3,1) S_true(3,3)], 'k', '--');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [S_mean(1,1) S_mean(1,3); S_mean(3,1) S_true(3,3)], 'k', '-');
hold on;
area.C95L = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [S_95L(1,1) S_95L(1,3); S_95L(3,1) S_95L(3,3)], 2/3*[1,1,1], '-');
hold on;
area.C95U = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [S_95U(1,1) S_95U(1,3); S_95U(3,1) S_95U(3,3)], 1/3*[1,1,1], '-.');
xlabel('X (mm)','fontsize', fsize);
ylabel('Z (mm)','fontsize', fsize);
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', fsize);
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeFLEXZPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeFLEXZPlot.jpg');

fprintf('FLE XZ Plane Area Perc. Diff.  Low: %f, Up: %f\n', ...
    100*(area.C95L-area.true)/area.true, 100*(area.C95U-area.true)/area.true);

figure(3);
clf;
% XY -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_true(1:2,1:2), 'k', '--');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_mean(1:2,1:2), 'k', '-');
hold on;
area.C95L = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_95L(1:2,1:2), 2/3*[1,1,1], '-');
hold on;
area.C95U = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', S_95U(1:2,1:2), 1/3*[1,1,1], '-.');
xlabel('X (mm)','fontsize', fsize);
ylabel('Y (mm)','fontsize', fsize);
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', fsize);
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeFLEXYPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeFLEXYPlot.jpg');

fprintf('FLE XY Plane Area Perc. Diff. Low: %f, Up: %f\n', ...
    100*(area.C95L-area.true)/area.true, 100*(area.C95U-area.true)/area.true);

figure(4);
clf;
% 3D view.
vol.true = plotCovarianceEllipse3RGB(conf, 50, [0 0 0]', S_true, 'k', 0);
%plotCovarianceEllipse3(conf, 50, [0 0 0]', S_mean, 'b', 1);
hold on;
vol.C95L = plotCovarianceEllipse3RGB(conf, 50, [-1 0 0]', S_95L, 2/3*[1,1,1], 1);
hold on;
vol.C95U = plotCovarianceEllipse3RGB(conf, 50, [1 0 0]', S_95U, 1/3*[1,1,1], 1);
xlabel('Y (mm)','fontsize', fsize);  
ylabel('Z (mm)','fontsize', fsize);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', fsize);  set(gca,'ZDir','reverse');
axis equal;
axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', fsize);
grid on;
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeFLE3DPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeFLE3DPlot.jpg');

fprintf('FLE Volume perc. diffs.  Low: %f, Up: %f\n', ...
    100*(vol.C95L-vol.true)/vol.true, 100*(vol.C95U-vol.true)/vol.true);


[TRE_true.RMS, TRE_true.Sigma, SigmaCovPA] = calcTRE(S_true, [mrk;tip0],0);
[TRE_mean.RMS, TRE_mean.Sigma, SigmaCovPA] = calcTRE(S_mean, [mrk;tip0],0);
[TRE_95L.RMS, TRE_95L.Sigma, SigmaCovPA] = calcTRE(S_95L, [mrk;tip0],0);
[TRE_95U.RMS, TRE_95U.Sigma, SigmaCovPA] = calcTRE(S_95U, [mrk;tip0],0);

figure(5);
clf;
% YZ -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_true.Sigma(2:3,2:3), 'k', '--');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_mean.Sigma(2:3,2:3), 'k', '-');
hold on;
area.C95L = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_95L.Sigma(2:3,2:3), 1/2*[1,1,1], '-');
hold on;
area.C95U = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_95U.Sigma(2:3,2:3), 1/3*[1,1,1], '-.');
xlabel('Y (mm)','fontsize', fsize);
ylabel('Z (mm)','fontsize', fsize);
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', fsize);
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeTREYZPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeTREYZPlot.jpg');

fprintf('TRE YZ Plane Area Perc. Diff.  Low: %f, Up: %f\n', ...
    100*(area.C95L-area.true)/area.true, 100*(area.C95U-area.true)/area.true);

figure(6);
clf;
% XZ -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [TRE_true.Sigma(1,1) TRE_true.Sigma(1,3); TRE_true.Sigma(3,1) TRE_true.Sigma(3,3)], 'k', '--');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [TRE_mean.Sigma(1,1) TRE_mean.Sigma(1,3); TRE_mean.Sigma(3,1) TRE_mean.Sigma(3,3)], 'k', '-');
hold on;
area.C95L = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [TRE_95L.Sigma(1,1) TRE_95L.Sigma(1,3); TRE_95L.Sigma(3,1) TRE_95L.Sigma(3,3)], 1/2*[1,1,1], '-');
hold on;
area.C95U = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', [TRE_95U.Sigma(1,1) TRE_95U.Sigma(1,3); TRE_95U.Sigma(3,1) TRE_95U.Sigma(3,3)], 1/3*[1,1,1], '-.');
xlabel('X (mm)','fontsize', fsize);
ylabel('Z (mm)','fontsize', fsize);
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', fsize);
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeTREXZPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeTREXZPlot.jpg');

fprintf('TRE XZ Plane Area Perc. Diff.  Low: %f, Up: %f\n', ...
    100*(area.C95L-area.true)/area.true, 100*(area.C95U-area.true)/area.true);

figure(7);
clf;
% XY -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_true.Sigma(1:2,1:2), 'k', '--');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_mean.Sigma(1:2,1:2), 'k', '-');
hold on;
area.C95L = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_95L.Sigma(1:2,1:2), 1/2*[1,1,1], '-');
hold on;
area.C95U = plotCovarianceEllipseRGB(conf, 1000/5, [0 0]', TRE_95U.Sigma(1:2,1:2), 1/3*[1,1,1], '-.');
xlabel('X (mm)','fontsize', fsize);
ylabel('Y (mm)','fontsize', fsize);
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', fsize);
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeTREXYPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeTREXYPlot.jpg');

fprintf('TRE YZ Plane Area Perc. Diff.  Low: %f, Up: %f\n', ...
    100*(area.C95L-area.true)/area.true, 100*(area.C95U-area.true)/area.true);

figure(8);
clf;
% 3D view.
vol.true = plotCovarianceEllipse3RGB(conf, 50, [0 0 0]', TRE_true.Sigma, 'k', 0);
%plotCovarianceEllipse3(conf, 50, [0 0 0]', TRE_mean.Sigma, 'b', 1);
hold on;
vol.C95L = plotCovarianceEllipse3RGB(conf, 50, [-1 0 0]', TRE_95L.Sigma, 2/3*[1,1,1], 1);
hold on;
vol.C95U = plotCovarianceEllipse3RGB(conf, 50, [1 0 0]', TRE_95U.Sigma, 1/3*[1,1,1], 1);
xlabel('Y (mm)','fontsize', fsize);  
ylabel('Z (mm)','fontsize', fsize);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', fsize);  set(gca,'ZDir','reverse');
axis equal;
axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', fsize);
grid on;
%set(gca, 'Visible', 'off');
print('-depsc2','-tiff','-r300', 'WestToolCovarianceRangeTRE3DPlot.eps');
print('-djpeg','-r300', 'WestToolCovarianceRangeTRE3DPlot.jpg');

fprintf('FLE Volume perc. diffs.  Low: %f, Up: %f\n', ...
    100*(vol.C95L-vol.true)/vol.true, 100*(vol.C95U-vol.true)/vol.true);
