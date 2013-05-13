function p=sdsvud(R,sig)
% sdsvud.m  2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% for 3-dim normal distribution with sigma(x)=sigma(y)=sigma(z)=sig
% find p=F(r)
% where F(r) = integral of f(x) from 0 to r.
% Maxwell pdf used for f(x) here.
%
if sig<=0 | R<0
  error(' check signs of input argument for cdcrud.m')
end

r=R/sig;
p=sf2(r);
