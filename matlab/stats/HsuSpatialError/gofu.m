%  gofu.m   09-26-97
%  g(u), Probability Density Function of the random variable u
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
function y=gofu(u)
global n

xm=sqrt(2*n-1);
fac1=1/(2^(n-1)*gamma(n/2));
y=fac1*(u+xm).^(n-1).*exp(-(u+xm).^2/4);
