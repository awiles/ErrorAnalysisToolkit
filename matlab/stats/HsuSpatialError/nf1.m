function y=nf1(x)
% nf1.m  3-06-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define f(x)=N(0,1)= Standard Normal Probability Density Function
% notice .^  is used  for vector x, sigma=1
fac1=1/sqrt(2*pi);
y=fac1*exp(-x.^2/2);
