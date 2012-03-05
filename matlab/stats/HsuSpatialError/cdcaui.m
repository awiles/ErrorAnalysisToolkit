function R=cdcaui(p,sig)
% cdcaui.m  1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given p and sigma, find R such that p = F(R),
% where F(R) = integral of f(x) from 0 to R.
%      Rayleigh pdf used for f(x) here.
%
if p>1 | p<0 | sig <=0
 error( 'check input arguments in cdcaui.m')
end

  if p==1
    R=Inf;
  else
    R=sig*sqrt(-2*log(1-p));
  end
