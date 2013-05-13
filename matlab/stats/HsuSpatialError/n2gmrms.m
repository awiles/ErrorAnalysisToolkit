function ratio=n2gmrms(n)
% n2gmrms.m  3-14-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
y=psin(n/2);
ratio=exp(y/2)/sqrt(n/2);
