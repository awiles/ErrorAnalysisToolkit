function z=bf1(r)
% bf1.m    r>=0,    1-15-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define f(r)= 2-dim normal pdf, rho=0, sigx and sigy may be unequal.
global mu

k1= r.^2*(1+mu^2)/(4*mu^2)        ;
k2= r.^2*(1-mu^2)/(4*mu^2)        ;
z= r/mu.* exp(-k1).*besseli(0,k2) ;

%
%  or use (bf1a.m)
%  function y=sigxy(x)
%  global mu xlam
%  fac1=1/(2*pi*mu)%;
%  gf=(cos(x).^2+sin(x).^2/mu^2)%;
%  y=fac1*(1-exp(-xlam^2.*gf/2))./gf %;
