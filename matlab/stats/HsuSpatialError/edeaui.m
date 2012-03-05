function [a,b]=edeaui(p,v,sx,sy)
% edeaui.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% elliptical normal pdf, x and y are uncorrelated
%  1/(2*pi*sx*sy) exp[(-x^2/sx^2 - y^2/sy^2)/2]
% elliptical area of integration
%    x^2/a^2 + y^2/b^2 = 1
%
%  v=a/b
%

if sx< 0 | sy< 0 | v <0 | (sx==0 & sy==0)
 error('check input arguments in edeaui.m')
end

b=1;
a=v*b;
p0=edeaud(a,b,sx,sy);

%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:2:5
  while (p0 > p)
    b=b-(0.1)^j;
    a=v*b;
    p0=edeaud(a,b,sx,sy);
  end

  while (p0 < p)
    b=b+(0.1)^j;
    a=v*b;
    p0=edeaud(a,b,sx,sy);
  end
end
