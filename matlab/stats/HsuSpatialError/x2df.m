function y=x2df(x)
% x2df.m  12-30-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
global n

if n<2
   fac1=1/(2^(n/2+1)*gamma(n/2+1));
   y=fac1*x.^(n/2).*exp(-x/2);
else
   fac1=1/(2^(n/2)*gamma(n/2));
   y=fac1*x.^(n/2-1).*exp(-x/2);
end
