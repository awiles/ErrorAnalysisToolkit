function p=edraud(a,b,sx,sy)
% edraud.m  2-07-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% The probability of a
% 2D joint elliptical normal pdf with sx, sy, rho=0
% over a rectangular region:  |x| < a,   |y| < b.
%

if sx<0 | sy<0 | a<0 | b<0 | (sx==0 & sy==0)
   error('check signs of input argument')
end

if sx==0
  p=nf2(b/sy);
end

if sy==0
  p=nf2(a/sx);
end

if sx > 0 & sy > 0
  p=nf2(a/sx)*nf2(b/sy);
end
