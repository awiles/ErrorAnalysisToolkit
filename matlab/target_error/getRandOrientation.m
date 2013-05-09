function [R, q] = getRandOrientation()

% generate random components of phi, theta
phi = 2*pi*rand(1,1);
theta = acos(1 - 2*rand(1,1));

% vector portion.
nx = sin(theta) * cos(phi);
ny = sin(theta) * sin(phi);
nz = cos(theta);

% generate random component about the vector.
thetaprime = 2*pi*rand(1,1);

q0 = cos(thetaprime/2);
qx = sin(thetaprime/2)*nx;
qy = sin(thetaprime/2)*ny;
qz = sin(thetaprime/2)*nz;

q = [q0, qx, qy, qz];

R = quat2rm(q);

