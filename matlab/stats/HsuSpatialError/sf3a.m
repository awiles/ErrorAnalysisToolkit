function y=sf3a(r)
% sf3a.m  2-19-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given p,  find G(r)=0
% define G(r)= F(r)-p={integral of trinormal pdf over [0,r]} - p
global p
% y=sf2(r) - p;
fac1=sqrt(2/pi);
y=fac1*(-r.*exp(-r.^2/2))+erf(r/sqrt(2)) - p ;
