clear all;
close all;

%% set up FLE components
SigmaIso = ((1/3)^2)/3 * eye(3);
SigmaAniso = diag([1/99 1/99 1/11]);


%% set up probe marker coordinates.  Rigid body parameters.
% define West tool design parameters.


testcase = 'C';

switch(testcase)
    case {'A'}
        casename = 'Des_d_Comp_85mmTip';
        probeDesign = 'd';
        A = 71;
        B = 54;
        rho = 85;
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
        rho = 800;
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

%% set up reference marker coordinates.  Rigid body parameters.
% define West tool design parameters.
r = 32;
% West design in Fig. 2(e).
ref.Rigid.mrk = [ -r/2 r/2 0;
    -r/2 -r/2 0;
    r/2 -r/2 0;
    r/2 r/2 0 ];
% tip is located at
ref.Rigid.tip = [ 0 0 0 ];

%% perform the hypothesis test M times.
N = 1000;

getRandOrientation = 0;

if(getRandOrientation)
    bValidTest = 0;
    while(~bValidTest)
        % rotate markers into measured space.
        pos = [0 0 -1000];
        screw.axis = getRandSphereSurfacePoints(1);
        screw.angle = 45 * randn(1);
        xfrm.R = screw2rot(screw);
        xfrm.pos = pos;
        [testResult, probeIso] = simTRE(SigmaIso, N, probeIso, xfrm);
        fprintf('Isotropic Error: %s\n', testResult);
        if( strcmp(testResult, 'Accept H0') | strcmp(testResult, 'Reject H0'))
            bValidTest = 1;
        end
        [testResult, probeAniso] = simTRE(SigmaAniso, N, probeAniso, xfrm);
        fprintf('Anisotropic Error: %s\n', testResult);
        if( (strcmp(testResult, 'Accept H0') | strcmp(testResult, 'Reject H0'))...
                && bValidTest == 1)
            bValidTest = 1;
        else
            bValidTest = 0;
        end
        fprintf('Anisotropic Error: %s\n', testResult);
    end
else % use this one.
    xfrm.pos = [0 0 -1000];
    rotX = 5; rotY = 0; rotZ = 30;
    xfrm.R = getRotMatrixd([rotX, rotY, rotZ]);
    [testResult, probeIso] = simTRE(SigmaIso, N, probeIso, xfrm);
    fprintf('Isotropic Error: %s\n', testResult);
    [testResult, probeAniso] = simTRE(SigmaAniso, N, probeAniso, xfrm);
    fprintf('Anisotropic Error: %s\n', testResult);
end

%% compare the accuracy of the two.

fprintf('Isotropic Error Covariance Matrices\n');
probeIso.Actual.stats.Sigma
probeIso.Meas.stats.Sigma

fprintf('Anisotropic Error Covariance Matrices\n');
probeAniso.Actual.stats.Sigma
probeAniso.Meas.stats.Sigma

%% plot stochastic cases and overlay the theoretical computation.
fig3 = figure('Name', 'Tip Error Distribution: Isotropic Error');
plotStochastic(probeIso.Meas.error', [0 0 0]', probeIso.Actual.stats.Sigma, [0 0 0]', probeIso.Meas.stats.Sigma );
fig4 = figure('Name', 'Tip Error Distribution: Anisotropic Error');
plotStochastic(probeAniso.Meas.error', [0 0 0]', probeAniso.Actual.stats.Sigma, [0 0 0]', probeAniso.Meas.stats.Sigma );

fig5 = figure('Name', 'Tip Error Covariance: Isotropic Error');
plotCovariance4View(N, [0 0 0]', probeIso.Actual.stats.Sigma, 'g', '-', 0)
plotCovariance4View(N, [0 0 0]', probeIso.Meas.stats.Sigma, 'k', '--', 1)

fig6 = figure('Name', 'Tip Error Covariance: Anisotropic Error');
plotCovariance4View(N, [0 0 0]', probeAniso.Actual.stats.Sigma, 'g', '-', 0)
plotCovariance4View(N, [0 0 0]', probeAniso.Meas.stats.Sigma, 'k', '--', 1)


%% compare Fitzpatrick formulation to ours.
fig7 = figure('Name', 'Comparison of Covariance Surfaces for Isotropic and Anisotropic Error');
plotCovarianceEllipse3(0.95, N/5, [0 0 0]', probeIso.Actual.stats.Sigma, 'b', 1);
plotCovarianceEllipse3(0.95, N/5, [0 0 0]', probeAniso.Actual.stats.Sigma, 'g', 1);
axis([-1.5, 1.5, -2.5, 2.5, -0.5, 0.5]);
%axis equal;
subtitle = sprintf('Case %s: RMS_{iso} = %.2f, RMS_{ani} = %.2f', testcase, probeIso.Actual.stats.RMS, probeAniso.Actual.stats.RMS);
figFontSize = 18;
%title({'\fontsize{figFontSize} Comparison of Isotropic and Anisotropic FLE'});
xlabel('Y','fontsize',figFontSize); ylabel('Z','fontsize',figFontSize); zlabel('X','fontsize',figFontSize);
set(gca,'fontsize',figFontSize);
title({'Comparison of Isotropic and Anisotropic FLE';subtitle}, 'fontsize', figFontSize);
print('-depsc','-tiff','-r300', casename);
colormap('gray');
%print -dpng -r300 iso_v_aniso;
%print -djpeg -r300 iso_v_aniso;
copyfile('*.eps', 'E:\docs\research\phd\publications\journal_papers\IEEE\anisotropic_TRE' );
fprintf('%s\n',subtitle);
% Principal Axes
fprintf('Isotropic Covariance Matrix in Principal Axes:\n');
probeIso.Actual.stats.SigmaPA
fprintf('Anisotropic Covariance Matrix in Principal Axes:\n');
probeAniso.Actual.stats.SigmaPA

fprintf('Isotropic Covariance Matrix in World Coordinates:\n');
probeIso.Actual.stats.Sigma
fprintf('Anisotropic Covariance Matrix in World Coordinates:\n');
probeAniso.Actual.stats.Sigma
