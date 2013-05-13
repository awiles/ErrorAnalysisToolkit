function xfrm = getRandomXfrm(s)
% returns the random transform corresponding to the fiducial tracking error
% covariance matrix s.  s is a 5x5 or 6x6 matrix.

if(size(s,1)< 5 || size(s,1) > 6)
    error('Invalid size of input matrix to getRandomXfrm given.');
end

[V,D] = eig(s);
x = (randn(1, size(s,1)))*(D^0.5)*V';
xfrm.pos = x(1:3);
if(size(s,1) == 5)
    xfrm.rot = [0 x(4:5) 0];
else
    xfrm.rot = [0 x(4:6)];
end
xfrm.rot(1) = sqrt(1 - xfrm.rot(2)^2 - xfrm.rot(3)^2 - xfrm.rot(4)^2);