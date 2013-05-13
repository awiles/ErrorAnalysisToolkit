function p=edeacd(a,b,sx,sy,rho)
% edeacd.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1996
%  Copyright 1996-1997 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% function p=edeacd(a,b,sx,sy,rho)
% edeac= elliptical normal pdf over elliptical area, x,y are correlated
% elliptical normal pdf
%  1/[2*pi*sx*sy*sqrt(1-rho^2)]*
%    exp{-[x^2/sx^2-2*rho*x*y/(sx*sy)+y^2/sy^2]/(2(1-rho^2)}
% elliptical area of integration
%    x^2/a^2 + y^2/b^2 = 1
%
% example
% a=0.6; b=1.2; sx=.4; sy=1; rho=0.3
% p=edeacd(0.6, 1.2, 0.4, 1, 0.3)
%

if abs(rho)>1 | sx<=0 | sy<=0 | a<=0 | b<=0
   error('check input arguments in edeacd(a,b,sx,sy,rho)')
end

   if abs(rho)==1
     u=sy/sx;
     st=sqrt((1+u^2))*sx;
     t0=sqrt(1+u^2)*a*b/sqrt(a^2*u^2+b^2);
     p=nf2(t0/st);
   end

   if rho==0
      p=edeaud(a,b,sx,sy)     ;
   end

   if abs(rho)<1 & abs(rho) >0
    sig1=sx/a;
    sig2=sy/b;

    A=1/sig1^2;
    B=-2*rho/(sig1*sig2);
    C=1/sig2^2;

    [a1, c1]=abc2ac(A,B,C);

    s1=1/sqrt(a1);
    s2=1/sqrt(c1);

    f0=sqrt(1-rho^2);

    t1=s1*f0;
    t2=s2*f0;

    p=edcaud(1,t1,t2) ;
   end
