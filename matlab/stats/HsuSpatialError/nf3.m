function r=nf3(p)
% nf3.m  3-07-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given p and f(x), find r such that G(r)= F(r)-p = 0,
% where F(r) = integral of f(x) from -r to r,
% using inverse erf(x) relation.
%   f(x)= N(0,1) used here.
%
if p>1 | p<0
  disp(' input error for nf3(p)')
  disp(' 0 <= p <= 1')
else
     if p==1
        r=Inf;
     else
        r=sqrt(2)*erfinv(p);
     end
end
