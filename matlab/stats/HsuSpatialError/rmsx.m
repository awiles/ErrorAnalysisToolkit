function r = rmsx(x)
% rmsx.m  8-06-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
%  find   sqrt{[x(1)^2+x(2)^2+ ... + x(N)^2]/N}

N=length(x);
r =sqrt(sum(x.^2)/N);
