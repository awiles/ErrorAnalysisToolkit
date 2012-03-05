function [a,b,c]=edbvui(p,k1,k2,sxyz)
% edbvui.m   2-28-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Given the probability  sx, sy, sz, k1=a/b, k2=b/c, p
% Find a, b
%
% 3D joint elliptical normal pdf with sx, sy, sz, rho=0
% over a rectangular region:  -a < x < a ,  -b < y < b, -c < z < c.
%

b=max(sxyz);
a=k1*b;
c=b/k2;
abc=[a b c];
p0=edbvud(abc,sxyz);

%%%%%%%%%  searching for a, b ,c %%%%%%%%%%%%%%%%
for j=1:1:5
  while (p0 > p)
    b=b-(0.1)^j ;
    a=k1*b       ;
    c=b/k2       ;
abc=[a b c];
p0=edbvud(abc,sxyz);
  end

  while (p0 < p)
    b=b+(0.1)^j ;
    a=k1*b       ;
    c=b/k2       ;
abc=[a b c];
p0=edbvud(abc,sxyz);
  end
end
