function y=x2disc2(x)
% x2disc2.m  1-16-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define G(a)=F(a)-p
% where F(a) = [integral of f(x) from 0 to a]
%     x™-distribution of degree n+2 is used for f(x) here.
%
global  p n
% chi-square pdf f(t,n) is singular at t=0 for n<2
% when integrating
% use F(x,n)=g(x,n)+F(x,n+2) relationship  for n<2
% g(x,n)= (x/2)^(n/2).*exp(-x/2)/gamma(n/2+1);
% F(x,n)=Integral of f(t,n) from 0 to x
% f(t,n)=chi-square pdf with degree n
% x2cdf.m  for F(x,n) and g(x,n)+F(x,n+2)

 y=x2cdf(x) - p;
