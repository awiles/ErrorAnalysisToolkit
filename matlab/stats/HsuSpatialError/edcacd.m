function p=edcacd(R,sx,sy,rho)
%  edcacd.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Elliptical PDF, correlated, Circular region

if R<0 | sx <0 | sy<0 | abs(rho)>1 | (sx==0 & sy==0)
   error('check input arguments in edcacd.m')
end

if abs(rho)==1
   st=sqrt(sx^2+sy^2);
   p=nf2(R/st);
end

if abs(rho)==0
   p=edcaud(R,sx,sy);
end

if abs(rho)<1 & abs(rho) >0
   A= 1/(1-rho^2)/sx^2  ;
   B= -2*rho/(1-rho^2)/(sx*sy) ;
   C= 1/(1-rho^2)/sy^2  ;
   [a1, c1]=abc2ac(A,B,C);
   s1=1/sqrt(a1);
   s2=1/sqrt(c1);
   p=edcaud(R,s1,s2);
end
