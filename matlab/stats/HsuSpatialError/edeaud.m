function p=edeaud(a,b,sx,sy)
% edeaud.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% edea= elliptical normal pdf over elliptical area,
% elliptical normal pdf, x and y are uncorrelated.
%  1/(2*pi*sx*sy) exp[(-x^2/sx^2 - y^2/sy^2)/2]
% elliptical area of integration
%    x^2/a^2 + y^2/b^2 = 1
%
% example
% a=0.6, b=1.2, sx=.4, sy=1.0.
% p=edeaud(0.6, 1.2, 0.4, 1.0)
%

if sx< 0 | sy< 0 | a<=0 | b<=0 | (sx==0 & sy==0)
 error('check input arguments in edeaud.m')
end

s1=sx/a;
s2=sy/b;
mu=min(s1,s2)/max(s1,s2) ;
p=edcaud(1,s1,s2);
