function plotTTEMonteCarlo(data, stats, casename, axisLim)

conf = 0.95;

%% local coordinate frame.
if(gcf == 1)
figure(1);
else
    figure(gcf+1);
end
clf;
plotCovariance4View(200, stats.simulated.local.mu', stats.simulated.local.cov, 'g', '-', 0);
plotCovariance4View(200, stats.theory.local.mu', stats.theory.local.cov, 'k', '-.', 1);
set(gca,'fontsize',18);
if(nargin > 2)
    filename = sprintf('%s-covariance-local', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

S1 = stats.simulated.local.cov;
S2 = stats.theory.local.cov;
figure(gcf+1);
clf;
% YZ -- top view.
hold on;
area.simulated = plotCovarianceEllipseRGB(conf, 1000/5, stats.simulated.local.mu(2:3)', S1(2:3,2:3), 2/3*ones(1,3), '-');
hold on;
area.theory = plotCovarianceEllipseRGB(conf, 1000/5, stats.theory.local.mu(2:3)', S2(2:3,2:3), 'k', '-.');
xlabel('Y (mm)','fontsize', 18);
ylabel('Z (mm)','fontsize', 18);
axis equal;
%axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', 18);
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-localYZ', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% XZ -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, stats.simulated.local.mu([1,3])', [S1(1,1) S1(1,3); S1(3,1) S1(3,3)], 2/3*ones(1,3), '-');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, stats.theory.local.mu([1,3])', [S2(1,1) S2(1,3); S2(3,1) S1(3,3)], 'k', '-.');
xlabel('X (mm)','fontsize', 18);
ylabel('Z (mm)','fontsize', 18);
axis equal;
%axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', 18);
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-localXZ', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% XY -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, stats.simulated.local.mu(1:2)', S1(1:2,1:2), 2/3*ones(1,3), '-');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, stats.theory.local.mu(1:2)', S2(1:2,1:2), 'k', '-.');
xlabel('X (mm)','fontsize', 18);
ylabel('Y (mm)','fontsize', 18);
axis equal;
%axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', 18);
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-localXY', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% 3D view.
vol.true = plotCovarianceEllipse3RGB(conf, 200, stats.simulated.local.mu', S1, 2/3*ones(1,3), 0);
hold on;
vol.C95L = plotCovarianceEllipse3RGB(conf, 200, stats.theory.local.mu', S2, 'k', 1);
xlabel('Y (mm)','fontsize', 18);  
ylabel('Z (mm)','fontsize', 18);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', 18);  set(gca,'ZDir','reverse');
axis equal;
%axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', 18);
grid on;
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-local3D', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

%% global coordinate frame.
figure(gcf+1);
clf;
plotCovariance4View(200, stats.simulated.global.mu', stats.simulated.global.cov, 'g', '-', 0);
plotCovariance4View(200, stats.theory.global.mu', stats.theory.global.cov, 'k', '--', 1);
set(gca,'fontsize',18);
if(nargin > 2)
    filename = sprintf('%s-covariance-global', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

S1 = stats.simulated.global.cov;
S2 = stats.theory.global.cov;
figure(gcf+1);
clf;
% YZ -- top view.
hold on;
area.simulated = plotCovarianceEllipseRGB(conf, 1000/5, stats.simulated.global.mu(2:3)', S1(2:3,2:3), 2/3*ones(1,3), '-');
hold on;
area.theory = plotCovarianceEllipseRGB(conf, 1000/5, stats.theory.global.mu(2:3)', S2(2:3,2:3), 'k', '-.');
xlabel('Y (mm)','fontsize', 18);
ylabel('Z (mm)','fontsize', 18);
axis equal;
%axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', 18);
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-globalYZ', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% XZ -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, stats.simulated.global.mu([1,3])', [S1(1,1) S1(1,3); S1(3,1) S1(3,3)], 2/3*ones(1,3), '-');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, stats.theory.global.mu([1,3])', [S2(1,1) S2(1,3); S2(3,1) S1(3,3)], 'k', '-.');
xlabel('X (mm)','fontsize', 18);
ylabel('Z (mm)','fontsize', 18);
axis equal;
%axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', 18);
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-globalXZ', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% XY -- top view.
hold on;
area.true = plotCovarianceEllipseRGB(conf, 1000/5, stats.simulated.global.mu(1:2)', S1(1:2,1:2), 2/3*ones(1,3), '-');
hold on;
area.mean = plotCovarianceEllipseRGB(conf, 1000/5, stats.theory.global.mu(1:2)', S2(1:2,1:2), 'k', '-.');
xlabel('X (mm)','fontsize', 18);
ylabel('Y (mm)','fontsize', 18);
axis equal;
%axis([-1.2 1.2 -1.2 1.2]);
view( [90 90] );
set(gca,'fontsize', 18);
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-globalXY', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% 3D view.
vol.true = plotCovarianceEllipse3RGB(conf, 200, stats.simulated.global.mu', S1, 2/3*ones(1,3), 0);
hold on;
vol.C95L = plotCovarianceEllipse3RGB(conf, 200, stats.theory.global.mu', S2, 'k', 1);
xlabel('Y (mm)','fontsize', 18);  
ylabel('Z (mm)','fontsize', 18);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', 18);  set(gca,'ZDir','reverse');
axis equal;
%axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', 18);
grid on;
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-global3D', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% 3D view.
plotCovarianceRings(conf, stats.simulated.global.mu', S1, 2/3*ones(1,3), '-', 2, [2 3 1]);
hold on;
plotCovarianceRings(conf, stats.theory.global.mu', S2, 'k', '-.', 2, [2 3 1]);
xlabel('Y (mm)','fontsize', 18);  
ylabel('Z (mm)','fontsize', 18);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', 18);  set(gca,'ZDir','reverse');
axis equal;
%axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', 18);
grid on;
%set(gca, 'Visible', 'off');
if(nargin > 2)
    filename = sprintf('%s-TTEcov-global3DRings', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end


figure(gcf+1);
clf;
% 3D view.
vol.C95L = plotCovarianceEllipse3RGB(conf, 200, stats.theory.global.mu', S2, 1/3*ones(1,3), 0);
xlabel('Y (mm)','fontsize', 28);  
ylabel('Z (mm)','fontsize', 28);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', 28);  set(gca,'ZDir','reverse');
axis equal;
%axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', 28);
grid on;
%set(gca, 'Visible', 'off');
if(nargin > 2)
    if(nargin > 3)
        axis(axisLim)
    end
    filename = sprintf('%s-TTEcov-global3DTheory', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

figure(gcf+1);
clf;
% 3D view.
plotCovarianceRings(conf, stats.theory.global.mu', S2, 'k', '-', 2, [2 3 1]);
xlabel('Y (mm)','fontsize', 28);  
ylabel('Z (mm)','fontsize', 28);  set(gca,'YDir','reverse');
zlabel('X (mm)','fontsize', 28);  set(gca,'ZDir','reverse');
axis equal;
%axis([-0.5 0.5 -1.2 1.2 -1.2 1.2]);
set(gca,'fontsize', 28);
grid on;
%set(gca, 'Visible', 'off');
if(nargin > 2)
    if(nargin > 3)
        axis(axisLim)
    end
    filename = sprintf('%s-TTEcov-global3DTheoryRings', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end


%% distribution.
figure(gcf+1);
[mu, sigma, rms, t] = plotProbHist(sqrt(sum(data.ttelocal.^2,2)), [], 'cumsum');
hold on;
plotDistance3DDist(stats.simulated.local.cov, t, '-', 2/3*ones(1,3));
plotDistance3DDist(stats.theory.local.cov, t, '-.', 'k');
hold off;
legend('TTE Distance Error Data', 'Simulated Covariance', 'Theoretical Covariance',...
    'Location', 'SouthEast');
xlabel('Distance Error (mm)', 'fontsize', 18);
ylabel('Probability', 'fontsize', 18);
title('Cummulative Distribution Histogram for TTE Distance Errors', 'fontsize', 18);
set(gca,'fontsize',18);

if(nargin > 2)
    filename = sprintf('%s-distance-distribution', casename);
    %print('-deps2','-tiff','-r300', [filename '.eps']);
    print('-dpng', '-r300', [filename '.png']);
end

max(t)