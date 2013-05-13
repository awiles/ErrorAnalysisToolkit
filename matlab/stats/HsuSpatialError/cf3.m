function r=cf3(p)
% cf3.m   3-06-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given p and f(x), find r such that G(r)=F(r)-p = 0,
% where F(r) = integral of f(x) from 0 to r.
%      Rayleigh pdf used for f(x) here.
%
if p>1 | p<0
  disp( 'input error for cf3(p)')
  disp( '0 <= p <= 1')
else
  if p==1
    r=Inf;
  else
    r=sqrt(-2*log(1-p));
  end
end
