function b=anl2b(a,N,lamb)
% anl2b.m   9-29-97
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given a, N and lambda,   find b

% if statistics toolbox is available
% ta=chi2inv(1-a,2*N-1);

% if statistics toolbox is not available, replace chi2inv with ichi2
  ta=ichi2(1-a,2*N-1);

Ka=sqrt(ta./(2*N))
tb=ta/lamb^2;

% if statistics toolbox is available
% b=chi2cdf(tb,2*N-1);

% if statistics toolbox is not available, replace chi2cdf with cdfchi2
  b=cdfchi2(tb,2*N-1);

