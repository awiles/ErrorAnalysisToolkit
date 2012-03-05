function R = quat2rm(q)
% insert reference here...
%#eml

w = 1; x = 2; y = 3; z = 4;
% compute the multiplications.
D = q' * q;

% build the rotation matrix.
R = zeros(3,3);
R(1,1) = D(w,w) + D(x,x) - D(y,y) - D(z,z);
R(1,2) = 2 * (-D(w,z) + D(x,y));
R(1,3) = 2 * ( D(w,y) + D(x,z));
R(2,1) = 2 * ( D(w,z) + D(x,y));
R(2,2) = D(w,w) - D(x,x) + D(y,y) - D(z,z);
R(2,3) = 2 * (-D(w,x) + D(y,z));
R(3,1) = 2 * (-D(w,y) + D(x,z));
R(3,2) = 2 * ( D(w,x) + D(y,z));
R(3,3) = D(w,w) - D(x,x) - D(y,y) + D(z,z);