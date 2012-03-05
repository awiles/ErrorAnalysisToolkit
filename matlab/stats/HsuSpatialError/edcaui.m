function R=edcaui(p,a,b)
% edcaui.m   1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Elliptical PDF, Uncorrelated, Circular region
% function R=edcaui(p,a,b)

if p<0 | p>1 | a< 0 | b< 0  | (a==0 & b==0)
    error('check signs of inputs for edcaui(p,a,b) ')
end

R=1;
p0=edcaud(R,a,b);
%%%%%%%%%%%%%%%%%%%%%%
for j=1:2:5
  while (p0 > p)
    R=R-(0.1)^j;
    p0=edcaud(R,a,b);
  end

  while (p0 < p)
    R=R+(0.1)^j;
    p0=edcaud(R,a,b);
  end
end
