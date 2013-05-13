function p=x2cdf(x)
% x2cdf.m  1-16-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define F(a)
% where F(a) = [integral of f(x) from 0 to a]
%     x™-distribution of degree n+2 is used for f(x) here.
%
global n
% chi-square pdf f(t,n) is singular at t=0 for n<2
% when integrating
% use F(x,n)=g(x,n)+F(x,n+2) relationship
% g(x,n)= (x/2)^(n/2).*exp(-x/2)/gamma(n/2+1);
%
% F(x,n)=Integral of f(t,n) from 0 to x
% f(t,n)=chi-square pdf with degree n
% x2df.m  defines f(t,n) and f(t,n+2)
% x2dg.m  defines g(x,n)

if n<2
  p=quad8('x2df',0,x) + x2dg(x,n);
else
  p=quad8('x2df',0,x);
end
