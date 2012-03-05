function a=cln2a(C,L,N)
% cln2a.m   9-29-97
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
     Ka=L/C;
     ta=2*N*Ka^2;
% if statistics toolbox is available
%    a=(1-chi2cdf(ta,2*N-1))*100 ;

% if statistics toolbox is not available, replace chi2cdf with cdfchi2
     a=(1-cdfchi2(ta,2*N-1))*100 ;
