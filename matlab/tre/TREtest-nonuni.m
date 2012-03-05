clear all;

%% set up the tool parameters.

% define the rigid body defintion in local probe coordinates.

%x = [0 0 0;
%    0 18.61 -36.95;
%   0 47.98 63.97];


% West design in Fig. 2(e).
A = 71;
B = 54;
rho = 170;

x = [ -A/2 B/2 0; ...
    -A/2 -B/2 0; ...
    A/2 -B/2 0; ...
    A/2 B/2 0 ];

% in drawing space, tip is located at
%tip0 = [-175 0 90];
tip0 = [ rho 0 0 ];

% now that the probe is created, let's rotate the probe into a useful
% orientation.
rotZ = 15;
R = getRotMatrixd([0, 0, rotZ]);
x = (R * x')';
tip0 = (R * tip0')';

% define an arbitrary FLE error value where each axis is identical.
FLE = 0.33;
TRE = calcTRE(FLE,[x;tip0]);
