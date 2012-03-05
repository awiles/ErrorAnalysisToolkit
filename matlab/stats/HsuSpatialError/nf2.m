function p=nf2(r)
% nf2.m  3-06-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given r>=0 , find p=F(r)
% F(r) = integral of N(0,1) over [-r ,r].   (using erf(x) relation)

if r<0
  error('check input argument in nf2(r) ')
end

p=erf(r/sqrt(2));
