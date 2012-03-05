function R=edcaci(p,sx,sy,rho)
% edcaci.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Elliptical PDF, correlated, Circular region

if p<0 | p>1 | sx <0 | sy<0 | abs(rho)>1 | (sx==0 & sy==0)
   error('check input arguments in edcacd.m')
end


R=1;
p0=edcacd(R,sx,sy,rho)
%%%%%%%%%%%%%%%%%%%%%%%
for j=1:2:5
  while (p0 > p)
    R=R-(0.1)^j;
    p0=edcacd(R,sx,sy,rho)
  end

  while (p0 < p)
    R=R+(0.1)^j;
    p0=edcacd(R,sx,sy,rho)
  end
end
