clear all;
close all;

cd e:\temp\TRE

%% set up FLE components
SigmaIso = ((1/3)^2)/3 * eye(3);
SigmaAniso = diag([1/99 1/99 1/11]);


%% set up probe marker coordinates.  Rigid body parameters.
% define West tool design parameters.


testcase = 'A';

switch(testcase)
    case {'A'}
        casename = 'Des_d_Comp_85mmTip';
        probeDesign = 'd';
        A = 71;
        B = 54;
        rho = 850;
        % define West tool design parameters.
        r = 32;
    case {'B'}
        casename = 'Des_d_Comp_170mmTip';
        probeDesign = 'd';
        Sigma = SigmaIso;
        A = 71;
        B = 54;
        rho = 170;
        % define West tool design parameters.
        r = 32;
    case {'C'}
        casename = 'Des_e_Comp_85mmTip';
        probeDesign = 'e';
        Sigma = SigmaIso;
        A = 71;
        B = 54;
        rho = 85;
        % define West tool design parameters.
        r = 32;
    case {'D'}
        casename = 'Des_e_Comp_170mmTip';
        probeDesign = 'e';
        Sigma = SigmaIso;
        A = 71;
        B = 54;
        rho = 170;
        % define West tool design parameters.
        r = 32;
    otherwise
        error('Ummm... that case does not exist');
end



switch(probeDesign)
    case {'d'}
        % West design in Fig. 2(d).
        probeIso.Rigid.mrk = [ 0 B/2 0;
            -A/2 0 0;
            0 -B/2 0;
            A/2 0 0 ];
        probeIso.Rigid.normals = repmat([ 0 0 1], 4, 1);
    case{'e'}
        % West design in Fig. 2(e).
        probeIso.Rigid.mrk = [ -A/2 B/2 0;
            -A/2 -B/2 0;
            A/2 -B/2 0;
            A/2 B/2 0 ];
        probeIso.Rigid.normals = repmat([ 0 0 1], 4, 1);
    otherwise
        error('Invalid tool design case.');
end
% tip is located at
probeIso.Rigid.tip = [ rho 0 0 ];



%copy to anisotropic data structure
probeAniso.Rigid.mrk = probeIso.Rigid.mrk;
probeAniso.Rigid.tip = probeIso.Rigid.tip;
probeAniso.Rigid.normals = probeIso.Rigid.normals;

% transform the bodies.
xfrm.pos = [0 0 -1000];
rotX = 0; rotY = 0; rotZ = 0;
xfrm.R = getRotMatrixd([rotX, rotY, rotZ]);
probeIso.Actual.mrk = (xfrm.R * probeIso.Rigid.mrk')' + repmat(xfrm.pos, 4, 1);
probeIso.Actual.normals = (xfrm.R * probeIso.Rigid.normals')';
probeIso.Actual.tip = (xfrm.R * probeIso.Rigid.tip')' + xfrm.pos;

probeAniso.Actual.mrk = (xfrm.R * probeAniso.Rigid.mrk')' + repmat(xfrm.pos, 4, 1);
probeAniso.Actual.normals = (xfrm.R * probeAniso.Rigid.normals')';
probeAniso.Actual.tip = (xfrm.R * probeAniso.Rigid.tip')' + xfrm.pos;

%% compute the TRE stats.
[probeIso.Actual.stats.RMS, probeIso.Actual.stats.Sigma, probeIso.Actual.stats.SigmaPA] =...
    calcTRE(SigmaIso, [probeIso.Actual.mrk; probeIso.Actual.tip], 0);
probeIso.Fitz.stats.RMS = calcTREold(sqrt(trace(SigmaIso)), [probeIso.Actual.mrk; probeIso.Actual.tip]);

[probeAniso.Actual.stats.RMS, probeAniso.Actual.stats.Sigma, probeAniso.Actual.stats.SigmaPA] =...
    calcTRE(SigmaAniso, [probeAniso.Actual.mrk; probeAniso.Actual.tip], 0);
probeAniso.Fitz.stats.RMS = calcTREold(sqrt(trace(SigmaAniso)), [probeAniso.Actual.mrk; probeAniso.Actual.tip]);

fprintf('Running Isotropic Simulation...\n');
[probeIso.Meas.error, probeIso.Meas.tip] = simTRE(SigmaIso, 10^5, probeIso.Rigid, probeIso.Actual );
probeIso = computeStats(probeIso);
fprintf('Running Anisotropic Simulation...\n');
[probeAniso.Meas.error, probeAniso.Meas.tip] = simTRE(SigmaAniso, 10^5, probeAniso.Rigid, probeAniso.Actual );
probeAniso = computeStats(probeAniso);


%% compare the accuracy of the two.
fprintf('Isotropic Error Covariance Matrices\n');
probeIso.Actual.stats.Sigma

fprintf('Anisotropic Error Covariance Matrices\n');
probeAniso.Actual.stats.Sigma

fig5 = figure('Name', 'Tip Error Covariance: Isotropic Error');
plotCovariance4View(1000, [0 0 0]', probeIso.Actual.stats.Sigma, 'g', '-', 0)
plotCovariance4View(1000, [0 0 0]', probeIso.Meas.stats.Sigma, 'k', '--', 1)

fig6 = figure('Name', 'Tip Error Covariance: Anisotropic Error');
plotCovariance4View(1000, [0 0 0]', probeAniso.Actual.stats.Sigma, 'g', '-', 0)
plotCovariance4View(1000, [0 0 0]', probeAniso.Meas.stats.Sigma, 'k', '--', 1)

%% compare Fitzpatrick formulation to ours.
fig7 = figure('Name', 'Comparison of Covariance Surfaces for Isotropic and Anisotropic Error');
plotCovarianceEllipse3(0.95, 1000/5, [0 0 0]', probeIso.Actual.stats.Sigma, 'b', 1);
plotCovarianceEllipse3(0.95, 1000/5, [0 0 0]', probeAniso.Actual.stats.Sigma, 'g', 1);
axis([-1.5, 1.5, -2.5, 2.5, -0.5, 0.5]);
%axis equal;
subtitle = sprintf('Case %s: RMS_{iso} = %.2f, RMS_{ani} = %.2f', testcase, probeIso.Actual.stats.RMS, probeAniso.Actual.stats.RMS);
figFontSize = 18;
%title({'\fontsize{figFontSize} Comparison of Isotropic and Anisotropic FLE'});
xlabel('Y','fontsize',figFontSize); ylabel('Z','fontsize',figFontSize); zlabel('X','fontsize',figFontSize);
set(gca,'fontsize',figFontSize);
title({'Comparison of Isotropic and Anisotropic FLE';subtitle}, 'fontsize', figFontSize);
%print('-depsc','-tiff','-r300', casename);
%print -dpng -r300 iso_v_aniso;
%print -djpeg -r300 iso_v_aniso;
%copyfile('*.eps', 'E:\docs\research\phd\publications\journal_papers\IEEE\anisotropic_TRE' );
fprintf('%s\n',subtitle);
% Principal Axes
fprintf('Isotropic Covariance Matrix in Principal Axes:\n');
probeIso.Actual.stats.SigmaPA
printforLATEX(probeIso.Actual.stats.SigmaPA, 'isoPA.txt')
fprintf('Anisotropic Covariance Matrix in Principal Axes:\n');
probeAniso.Actual.stats.SigmaPA
printforLATEX(probeAniso.Actual.stats.SigmaPA, 'anisoPA.txt')

fprintf('Isotropic Covariance Matrix in World Coordinates:\n');
probeIso.Actual.stats.Sigma
printforLATEX(probeIso.Actual.stats.Sigma, 'iso.txt')
fprintf('Anisotropic Covariance Matrix in World Coordinates:\n');
probeAniso.Actual.stats.Sigma
printforLATEX(probeAniso.Actual.stats.Sigma, 'aniso.txt')
