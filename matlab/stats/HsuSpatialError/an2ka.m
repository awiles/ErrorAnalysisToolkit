function Ka=an2ka(a,N)
% an2ka.m  9-29-97
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given a and N,  find Ka

% if statistics toolbox is available, use the command chi2inv.m
%   ta=chi2inv(1-a,2*N-1);

% if statistics toolbox is not available, to replace chi2inv.m
    ta=ichi2(1-a,2*N-1);

Ka=sqrt(ta/(2*N));
