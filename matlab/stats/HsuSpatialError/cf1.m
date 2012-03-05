function y=cf1(r)
% cf1.m  6/04/97,   define f(r)= circular normal pdf
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define f(r)= circular normal pdf
% Rayleigh Probability Density Function  sigma=1
y=r.* exp(-r.^2/2);
