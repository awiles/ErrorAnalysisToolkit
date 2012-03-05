function [p_close,dist] = getPointtoPlane(n, orig, pt) 

%normalize the normal.
n = 1/(sum(n.^2)^(1/2)) * n;

dist = n * (pt-orig)';

p_close = pt - dist*n;
