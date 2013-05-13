function p=cdcaud(R,sig)
% cdcaud.m  1-23-96
%
%  Spatial Error Analysis ToolBoxc Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% r=R/sig
% for 2-dim normal distribution with sigma(x)=sigma(y)=sig
% find p=F(r)
% where F(r) = integral of f(x) from 0 to r.
% Rayleigh pdf used for f(x) here.
%
if sig<=0 | R<0
  error(' check signs of input argument for cdcaud.m')
end

r=R/sig;
  p=1-exp(-r.^2/2);
end
