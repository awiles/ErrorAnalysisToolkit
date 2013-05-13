function p=gf2(r,mu)
% gf2.m  01-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% mu=min(sx,sy)/max(sx,sy)
% r=R/max(sx,sy)
%

% r>=0, mu in (0,1)
 if r<0 | mu<=0 | mu>=1
   error(' check input arguments in gf2(r,mu)')
 end

   n=2000;
   h=2*pi/n;
   fac1=1/(2*pi*mu);
   x=[0:n]*h;
   gf=(cos(x).^2+sin(x).^2/mu^2);
   y=(1-exp(-r^2.*gf/2))./gf;
   p=fac1*simprule(y,h);
