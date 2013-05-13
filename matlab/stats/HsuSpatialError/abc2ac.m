function [a1,c1]=abc2ac(a,b,c)
% abc2ac.m  7-10-97
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% rotation of ellipse
% Ax^2 + Bxy + Cy^2 = 1   ---> au^2 + cv^2 = 1.

ang=atan2(b,a-c);
t=ang/2;
tdeg=t*180/pi;

ct=cos(t);
st=sin(t);

a1=a*ct^2+b*ct*st+c*st^2;
c1=a*st^2-b*ct*st+c*ct^2;
