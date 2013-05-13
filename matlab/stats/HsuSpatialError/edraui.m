function [a,b]=edraui(p,v,sx,sy)
% edraui.m   2-07-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Given the probability  sx, sy, v=a/b, p
% Find a, b
%
% 2D joint elliptical normal pdf with sx, sy, rho=0
% over a rectangular region:  -a < x < a ,  -b < y < b.
%

b=max(sx,sy);
a=v*b;
p0=edraud(a,b,sx,sy);

%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:2:5
  while (p0 > p)
    b=b-(0.1)^j ;
    a=v*b       ;
    p0=edraud(a,b,sx,sy);
  end

  while (p0 < p)
    b=b+(0.1)^j ;
    a=v*b       ;
    p0=edraud(a,b,sx,sy);
  end
end
