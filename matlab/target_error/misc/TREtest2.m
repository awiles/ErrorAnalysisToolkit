% TREtest2.m
close all;
clear all;

% polaris active case 3 on pg. 541(9) of West.
A = 71;
B = 54;
rho = 170;
r = 64;
d = 200;

% probe marker and tip positions
xp = [-A/2 -B/2 0; -A/2 B/2 0; A/2 B/2 0; A/2 -B/2 0];
tp = [rho 0 0];

% coordinate reference frame
ang = deg2rad(45);
xc = [r*cos(ang) r*sin(ang) 0; -r*cos(ang) r*sin(ang) 0;...
        -r*cos(ang) -r*sin(ang) 0; r*cos(ang) -r*sin(ang) 0];
tc = [d 0 0];

% define an arbitrary FLE error value where each axis is identical.
FLE = 0.33;
TREp_West = calcTRE_NDI(FLE,[xp;tp],1);
TREc_West = calcTRE_NDI(FLE,[xc;tc],1);
TREpc_West = sqrt(TREp_West^2 + TREc_West^2);

% redefine the FLE error to have 3 times the error on the z-axis.
FLE = [0.0995 0.0995 0.2985].^2;
TRE0 = calcTRE_NDI(FLE,[xp;tp],2);

ang = [-45:5:45];

for i = 1:length(ang)
    R = getRotMatrix(0, ang(i), 0)
    xp_test = (R * xp')';
    tp_test = (R * tp')';
    xc_test = (R * xc')';
    tc_test = (R * tc')';
    
    TREp_ang(i) = calcTRE_NDI(FLE,[xp_test;tp_test],0);
    TREc_ang(i) = calcTRE_NDI(FLE,[xc_test;tc_test],0);
    TREpc_ang(i) = sqrt(TREp_ang(i)^2 + TREc_ang(i));
end

figure(3);
plot(ang, TREp_ang);
hold on;
plot(ang, TREc_ang,'g');
plot(ang, TREpc_ang,'k');
plot(0, TREp_West, 'ob');
plot(0, TREc_West, 'og');
plot(0, TREpc_West, 'ok');
hold off;
xlabel('Rotation about y-axis from default orientation (degrees)');
ylabel('TRE (mm)');
legend('TRE - Probe', 'TRE - CRF', 'TRE - Combined',0);

fprintf('Probe TRE = %4.2f\nCRF TRE = %4.2f\nCombined TRE = %4.2f\n', TREp_West, TREc_West, TREpc_West);