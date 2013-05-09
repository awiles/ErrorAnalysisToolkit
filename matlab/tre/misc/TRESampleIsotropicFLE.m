% TREscript.m
%**************************************************************************
% Target Registration Error Testing Script
% 
%       Written by Andrew Wiles, March 24, 2005
%       Copyright (C)2005, Northern Digital Inc.  All rights reserved.
%
%       Confidential and Proprietary.
%       Not for distribution outside of NDI without proper authorization.
%
%**************************************************************************
%   This script is used to calcuate the TRE for a given tool definition.
%   Three simulations are completed here:
%       Test 1: Calculate the TRE for tip positions located along a line.
%       Test 2: Calculate the TRE for tip positions located along a circle
%       rotated about a given axis.
%       Test 3: Calculate the TRE for a given tip direction but different
%       lengths.
%**************************************************************************
clear all;

% define an arbitrary FLE error value.
FLE = 0.35;

% define the rigid body defintion.
x = [0 0 0;
    0 18.61 -36.95;
    0 47.98 63.97];

%**************************************************************************
% Test #1
%**************************************************************************
% in global space, tip is located at 
tip0 = [-175 0 90]';
% rotate into local space.
pitch = 50;
R = [cos(pitch) 0 -sin(pitch);
    0 1 0;
    sin(pitch) 0 cos(pitch)];
tip0 = R * tip0;

% determine TRE for tip at tip0 and varying in the z direction
zpos = [-100:5:100]';

% FIRST, calculate for local z-axis.
tip = repmat(tip0',size(zpos,1),1) + [zeros(size(zpos)) zeros(size(zpos)) zpos];
for i = 1:size(tip,1)
    TRE_z(i) = calcTRE(FLE,[x;tip(i,:)]);
end
TRE_z0 = calcTRE(FLE,[x;tip0']);

% SECOND, calculate for global z-axis.
tipG = repmat(tip0',size(zpos,1),1) + (R'*[zeros(size(zpos)) zeros(size(zpos)) zpos]')';
for i = 1:size(tip,1)
    TRE_zG(i) = calcTRE(FLE,[x;tipG(i,:)]);
end
TRE_z0 = calcTRE(FLE,[x;tip0']);

figure(1);
plot(tip(:,3),TRE_z);
hold on;
plot(tip(:,3),TRE_zG,'g');
scatter(tip0(3),TRE_z0,'o');
hold off;
title('TRE vs. tip position in Z-axis');
xlabel('Z axis position (mm)');
ylabel('TRE (mm)');
legend('TRE local Z-axis','TRE global Z-axis','Specific Case in Drawing',0);

%**************************************************************************
% Test #2
%**************************************************************************
% determine TRE for tips located at constant distance from the origin
% rotated about the y-axis.

radius = 100;                         %distance from origin                  
alpha = (pi/180)*([0:1:360]');     %angle about origin

% define tip positions on the unit circle and then set the radius.
tip_xy = radius * [ cos(alpha) sin(alpha) zeros(size(alpha)) ];
tip_xz = radius * [ cos(alpha) zeros(size(alpha)) sin(alpha) ];
tip_yz = radius * [ zeros(size(alpha)) cos(alpha) sin(alpha) ];

% cycle through and determine TRE for each case.
for i = 1:size(tip_xy,1)
    TRE_xy(i) = calcTRE(FLE,[x;tip_xy(i,:)]);
    TRE_xz(i) = calcTRE(FLE,[x;tip_xz(i,:)]);
    TRE_yz(i) = calcTRE(FLE,[x;tip_yz(i,:)]);
end

% plot results.
figure(2);
subplot(3,1,1);
polar(alpha, TRE_xy');
title('TRE as a function of angle about the Z-axis, R_{tip}=100mm');
subplot(3,1,2);
polar(alpha, TRE_xz');
title('TRE as a function of angle about the Y-axis, R_{tip}=100mm');
subplot(3,1,3);
polar(alpha, TRE_yz');
title('TRE as a function of angle about the X-axis, R_{tip}=100mm');

%**************************************************************************
% Test #3
%**************************************************************************

% determine the direction vector based on the arbitrary tip position given
% above, i.e. the drill guide tip position.
dirVec = 1/norm(tip0)*tip0;

% set which distances we are interested in.
distance = [0:5:300]';

% calculate the various tip positions in the direction of interest.
tip_dist = [(distance .* repmat(dirVec(1), size(distance,1),1)),...
        (distance .* repmat(dirVec(2), size(distance,1),1)),...
        (distance .* repmat(dirVec(3), size(distance,1),1))];

% cycle through and determine TRE for each case.
for i = 1:size(tip_dist,1)
    TRE_dist(i) = calcTRE(FLE,[x;tip_dist(i,:)]);    
end

% plot.
figure(3);
plot(distance,TRE_dist);
title('TRE vs. Tip Distance');
xstring = sprintf('Distance from origin in direction [%3.2f,%3.2f,%3.2f]', dirVec);
xlabel(xstring);
ylabel('TRE (mm)');