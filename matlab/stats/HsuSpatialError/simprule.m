function s=simprule(y,h)
% simprule.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
n= length(y)-1   ;
s0=y(1)+y(n+1);
s1=sum(y(2:2:n));
s2=sum(y(3:2:n-1));
s=h*(s0+4*s1+2*s2)/3;
