function p=edca(r,mu)
%  edca.m   1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% p  = probability of ED (elliptical normal pdf) over CA(circular area)
% mu = min(sx,sy)/max(sx,sy)
% r  = R/max(sx,sy)
% x  = theta angle
%      Normalized version of edcaud.m  { p=edcaud(R,sx,sy) }
% 1-23-96  from bf1b.m

     if mu>1 |  mu<0 | r < 0
     error('check input arguments in edca(r,mu) ')
     end

mlim=0.01;
if mu<mlim
% treated as one-dimensional case (mu==0)
  p=nf2(r);
end

if mu>=mlim  & mu<1
   p=gf2(r,mu);   % gf2.m  replaces the next 7 lines
end

if mu==1
  p=cf2(r);
end

% method 1
% this part does not work for mu=.5, and r=4
%   clear global
%   global mu r
%   p=quad8('bf1a',0,2*pi);
%   insure
% this part does not work for mu=.5, and r=4
%
% method 2
% this part works nicely, and faster
%  p=gf2(r,mu)    % gf2.m  replaces the next 7 lines
%    n=2000;
%    h=2*pi/n;
%    fac1=1/(2*pi*mu);
%    x=[0:n]*h;
%    gf=(cos(x).^2+sin(x).^2/mu^2);
%    y=(1-exp(-r^2.*gf/2))./gf;
%    p=fac1*simprule(y,h);
% this part works nicely
