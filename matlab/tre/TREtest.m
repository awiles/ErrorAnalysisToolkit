clear all;

% define West tool design parameters.
A = 71;
B = 54;
rho = 170;

% define the rigid body defintion in local probe coordinates.

%x = [0 0 0;
%    0 18.61 -36.95;
%   0 47.98 63.97];

% West design in Fig. 2(e).
x = [ -A/2 B/2 0;
    -A/2 -B/2 0;
    A/2 -B/2 0;
    A/2 B/2 0 ];


% in drawing space, tip is located at 
%tip0 = [-175 0 90];
tip0 = [ rho 0 0 ];

% now that the probe is created, let's rotate the probe into a useful
% orientation.
rotZ = 0;
R = getRotMatrixd([0, 0, rotZ]);
x = (R * x')';
tip0 = (R * tip0')';

% define an arbitrary FLE error value where each axis is identical.
FLE = 0.33;
TRE_West = calcTRE_NDI(FLE,[x;tip0],1);

% redefine the FLE error to have 3 times the error on the z-axis.
FLE = 0.33;
TRE0 = calcTRE_NDI(FLE,[x;tip0],1);

ang = [-45:1:45];

for i = 1:length(ang)
    R = getRotMatrixd([0, ang(i), 0]);
    x_test = (R * x')';
    tip0_test = (R * tip0')';
    
    TRE_ang(i) = calcTRE_NDI(FLE,[x_test;tip0_test]);
    TRE_ang0(i) = calcTRE(FLE, [x_test;tip0_test] );
end

FLE = [0.0995 0.0995 0.2985].^2;
TRE0 = calcTRE_NDI(FLE,[x;tip0],1);

ang = [-45:1:45];

for i = 1:length(ang)
    R = getRotMatrixd([0, ang(i), 0]);
    x_test = (R * x')';
    tip0_test = (R * tip0')';
    
    TRE_ang2(i) = calcTRE_NDI(FLE,[x_test;tip0_test]);
end

figure(3);
plot(ang, TRE_ang);
hold on;
plot(ang, TRE_ang2, 'g');
hold off;
xlabel('Rotation about y-axis from default orientation (degrees)');
ylabel('TRE (mm)');
legend('Isotropic \Sigma', 'Anisotropic \Sigma' );
title('Comparison of TRE Model using Isotropic and Anisotropic FLE');
