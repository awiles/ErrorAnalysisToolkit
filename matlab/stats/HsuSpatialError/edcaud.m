function p=edcaud(R,sx,sy)
% edcaud.m 1-23-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Elliptical PDF, Uncorrelated, Circular region
% edcaud.m, Find probability  p, given R,sx,sy.
% bf1b.m is the normalized version
% function p=edcaud(R,sx,sy)
%   1>= mu >=0

if R<0 | sx< 0 | sy< 0 | (sx==0 & sy==0)
    error('check input  arguments for edcaud(R,sx,sy) ')
end

mu=min(sx,sy)/max(sx,sy);
r=R/max(sx,sy);
p=edca(r,mu);
