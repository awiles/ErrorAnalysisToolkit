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

% random fiducials
probe = generateRandFiducialTargets(3, 100, 200);
x = probe.Rigid.mrk
tip0 = probe.Rigid.tip

% x = [26.1085   34.4732  -10.4158; ...
%    29.8527  -44.1469  -35.7065; ...
%     1.4457   -8.9735  -35.2046 ];
% tip0 = [  -68.2192 39.5070    4.6551];

% now that the probe is created, let's rotate the probe into a useful
% orientation.
rotZ = 15*randn(1);
R = getRotMatrixd([0, 0, rotZ]);
x = (R * x')';
tip0 = (R * tip0')';

% define an arbitrary FLE error value where each axis is identical.
FLE1 = 0.33;
[TRE(1) cov covPA] = calcTRE(FLE1,[x;tip0], 0);
fprintf('Test 1: TRE = %f\n', TRE(1));
cov

FLE2 = (0.33)^2/3*eye(3);
[TRE(2) cov covPA] = calcTRE(FLE2,[x;tip0], 0);
fprintf('Test 2: TRE = %f\n', TRE(2));
cov

FLE0 = (0.33)^2/3*eye(3);
FLE3 = zeros(3,3,size(x,1));
for i = 1:size(x,1)
    FLE3(:,:,i) = FLE0;
end
FLE3;
[TRE(3) cov covPA] = calcTRE(FLE3,[x;tip0], 0);
fprintf('Test 3: TRE = %f\n', TRE(3));
cov


FLE0 = (0.33)^2/3*eye(3);
FLE4 = zeros(3,3,size(x,1));
for i = 1:size(x,1)
    FLE4(:,:,i) = 0.5 *i * FLE0;
end
[TRE(4) cov covPA] = calcTRE(FLE4,[x;tip0], 0);
test.mrk = x;
test.tip = tip0;
[probe.Meas.error, probe.Meas.tip] = simTRE(FLE4, 100000, test, test );
probe = computeStats(probe);
fprintf('Test 4: TRE_t = %f, TRE_m = %f\n', TRE(4), probe.Meas.stats.RMS);
fprintf('Test #4 - theory\n');
cov
fprintf('Test #4 - simulation\n');
probe.Meas.stats.Sigma

FLE5 = weightMatrix(0.33^2, [1 1 3]);
[TRE(5) cov covPA] = calcTRE(FLE5,[x;tip0], 0);
fprintf('Test 5: TRE = %f\n', TRE(5));
cov

FLE6 = zeros(3,3,size(x,1));
for i = 1:size(x,1)
    FLE6(:,:,i) = weightMatrix(0.33^2, [1 1 3]);
end
[TRE(6) cov covPA] = calcTRE(FLE6,[x;tip0], 0);
fprintf('Test 6: TRE = %f\n', TRE(6));
cov

FLE7 = mean(FLE6,3);
[TRE(7) cov covPA] = calcTRE(FLE7,[x;tip0], 0);
fprintf('Test 7: TRE = %f\n', TRE(7));
cov
