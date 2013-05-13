function y=psi(x)
% psi.m  6-21-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
global n

if x==1
y=n;
else
y=(1-x.^n)./(1-x); % use 1-eps to avoid singularity when integrating
end
