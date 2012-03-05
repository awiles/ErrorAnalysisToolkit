function y=x2dg(x,n)
% x2dn.m  1-2-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
if n<2
fac1=1/(2^(n/2)*gamma(n/2+1));
y=fac1*x.^(n/2).*exp(-x/2);
else
 y=0;
end
