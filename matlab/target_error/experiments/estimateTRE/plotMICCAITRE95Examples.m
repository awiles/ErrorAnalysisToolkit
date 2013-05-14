%% script to plot the MICCAI 95 percent confidence intervals.

cd 'E:\docs\research\phd\publications\posters\MICCAI2007'

design = 'e';
A = 71;
B = 54;
rho = 85;
r = 64;
model = 'Aniso';

SigmaIso = ((1/3)^2)/3 * eye(3);
SigmaAniso = diag([1/99, 1/99, 1/11]);

if(strcmp(model, 'Iso'))
    Sigma = SigmaIso;
elseif(strcmp(model, 'Aniso'))
    Sigma = SigmaAniso;
end

filenameprefix = sprintf('%s_%03d', model, r);

[probe.Rigid.mrk, probe.Rigid.normals, probe.Rigid.tip] = getWestToolDesign(design, A, B, rho);
[ref.Rigid.mrk, ref.Rigid.normals, ref.Rigid.tip] = getWestToolDesign('e', r, r, 0);

%% set up the probe transform.
xfrm.pos = [100 100 -1000];
xfrm.R = getRotMatrixd([-10 -15 0]);
probe.xfrm = xfrm;
% transform rigid body into test space.
probe.Actual.mrk = (xfrm.R * probe.Rigid.mrk')' + repmat(xfrm.pos, 4, 1);
probe.Actual.normals = (xfrm.R * probe.Rigid.normals')';
probe.Actual.tip = (xfrm.R * probe.Rigid.tip')' + xfrm.pos;

%% set up the reference transform.
xfrm.pos = [-0 50 -950];
xfrm.R = getRotMatrixd([0 15 0]);
ref.xfrm = xfrm;
% transform rigid body into test space.
ref.Actual.mrk = (ref.xfrm.R * ref.Rigid.mrk')'...
    + repmat(ref.xfrm.pos, size(ref.Rigid.mrk,1), 1);
ref.Actual.normals = (ref.xfrm.R * ref.Rigid.normals')';
ref.Actual.tip = probe.Actual.tip;
ref.Rigid.tip = (ref.xfrm.R' * (ref.Actual.tip - ref.xfrm.pos)')';

ref_tip_dist = sqrt(sum((ref.Rigid.tip) .^2))

% compute errors.
probe.Actual.stats.mu = zeros(1,3);
[probe.Actual.stats.RMS, probe.Actual.stats.Sigma, probe.Actual.stats.SigmaPA] =...
    calcTRE(Sigma, [probe.Actual.mrk; probe.Actual.tip], 1);

ref.Actual.stats.mu = zeros(1,3);
[ref.Actual.stats.RMS, ref.Actual.stats.Sigma, ref.Actual.stats.SigmaPA] =...
    calcTRE(Sigma, [ref.Actual.mrk; ref.Actual.tip], 2);

combined.Actual.stats.mu = probe.Actual.stats.mu + ref.Actual.stats.mu;
combined.Actual.stats.RMS = sqrt(probe.Actual.stats.RMS^2 + ref.Actual.stats.RMS^2);
combined.Actual.stats.Sigma = ref.Actual.stats.Sigma + probe.Actual.stats.Sigma;
combined.Actual.stats.Sigma

% plot.
axis_bounds = [-1.5 1.5 -5.0 5.0 -0.5 0.5];
view_angles = [37.5 25];
figFontSize = 24;

figure(3); clf;
plotCovarianceEllipse3(0.95, 200, [0 0 0]', probe.Actual.stats.Sigma, 'b', 0);
axis(axis_bounds);
view(view_angles);
titlestring = sprintf('RMS = %2.2f', probe.Actual.stats.RMS);
title(titlestring, 'fontsize', figFontSize);
set(gca,'fontsize',figFontSize);
xlabel('Y','fontsize',figFontSize);
ylabel('Z', 'fontsize',figFontSize);
zlabel('X', 'fontsize',figFontSize);
filenameout = sprintf('%s_probe_only', filenameprefix);
print('-dpng', '-r300', filenameout);

figure(4); clf;
plotCovarianceEllipse3(0.95, 200, [0 0 0]', ref.Actual.stats.Sigma, 'r', 0);
axis(axis_bounds);
view(view_angles);
titlestring = sprintf('RMS = %2.2f', ref.Actual.stats.RMS);
title(titlestring, 'fontsize', figFontSize);
set(gca,'fontsize',figFontSize);
xlabel('Y','fontsize',figFontSize);
ylabel('Z', 'fontsize',figFontSize);
zlabel('X', 'fontsize',figFontSize);
filenameout = sprintf('%s_ref_only', filenameprefix);
print('-dpng', '-r300', filenameout);

figure(5); clf;
plotCovarianceEllipse3(0.95, 200, [0 0 0]', combined.Actual.stats.Sigma, 'g', 0);
axis(axis_bounds);
view(view_angles);
titlestring = sprintf('RMS_{combined} = %2.2f', combined.Actual.stats.RMS);
title(titlestring, 'fontsize', figFontSize);
set(gca,'fontsize',figFontSize);
xlabel('Y','fontsize',figFontSize);
ylabel('Z', 'fontsize',figFontSize);
zlabel('X', 'fontsize',figFontSize);
filenameout = sprintf('%s_combined_probes', filenameprefix);
print('-dpng', '-r300', filenameout);

figure(6); clf;
plotCovarianceEllipse3(0.95, 200, [0 0 0]', probe.Actual.stats.Sigma, 'b', 0);
plotCovarianceEllipse3(0.95, 200, [0 0 0]', ref.Actual.stats.Sigma, 'r', 1);
axis(axis_bounds);
view(view_angles);
titlestring = sprintf('RMS_{probe} = %2.2f, RMS_{ref} = %2.2f', probe.Actual.stats.RMS, ref.Actual.stats.RMS);
title(titlestring, 'fontsize', figFontSize);
set(gca,'fontsize',figFontSize);
xlabel('Y','fontsize',figFontSize);
ylabel('Z', 'fontsize',figFontSize);
zlabel('X', 'fontsize',figFontSize);
%subtitle = sprintf('Isotropic FLE: RMS_{probe} = %2.2f, RMS_{ref} = %2.2f',...
%    probe.Actual.stats.RMS, ref.Actual.stats.RMS);
%title({'Probe and Reference Overlayed';subtitle});
filenameout = sprintf('%s_ref+probe', filenameprefix);
print('-dpng', '-r300', filenameout);