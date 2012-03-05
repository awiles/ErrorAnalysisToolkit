function p=cf2(r)
% cf2.m   1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% for 2-dim normal distribution with sigma(x)=sigma(y)
% find p=F(r)
% where F(r) = integral of f(x) from 0 to r.
% Rayleigh pdf used for f(x) here.
%

if r<0
  error('check input argument in cf2(r)')
end

p=1-exp(-r.^2/2);
