function y=tdis(x)
%  tdis.m    3-22-94
%
%  Spatial Error Analysis ToolBoxc Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% student t - Probability Density Function
% notice .* and .^  are used  for vector x:
global n
fac1=gamma((n+1)/2)/gamma(n/2)/sqrt(pi*n) ;
y=fac1*(1+x.^2/n).^(-(n+1)/2);
