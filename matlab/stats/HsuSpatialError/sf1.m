function y=sf1(r)
% sf1.m    3-06-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define f(r)= Maxwell Probability Density Function  sigma=1
y=sqrt(2/pi)*r.^2 .* exp(-r.^2/2);
