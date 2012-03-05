function [a,b]=edeaci(p,v,sx,sy,rho)
% edeaci.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% function [a,b]=edeaci(p,v,sx,sy,rho)
% edeac= elliptical normal pdf over elliptical area, x,y are correlated
% elliptical normal pdf
%  1/[2*pi*sx*sy*sqrt(1-rho^2)]
%     exp{-[x^2/sx^2-2*rho*x*y/(sx*sy)+y^2/sy^2]/(2(1-rho^2)}
% elliptical area of integration
%    x^2/a^2 + y^2/b^2 = 1
%
%           v=a/b
%
% example
% p=0.7; v=0.5; sx=.4; sy=1; rho=0.3
% [a,b]=edeaci(0.7, 0.5, 0.4, 1, 0.3)
%

if abs(rho)>1 | v<=0 |  sx<=0 | sy<=0
   error('check input arguments in edeacd.m ')
end


b=1;
a=v*b;
p0=edeacd(a,b,sx,sy,rho);

%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:2:5
  while (p0 > p)
    b=b-(0.1)^j;
    a=v*b;
    p0=edeacd(a,b,sx,sy,rho);
  end

  while (p0 < p)
    b=b+(0.1)^j;
    a=v*b;
    p0=edeacd(a,b,sx,sy,rho);
  end
end
