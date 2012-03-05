function y=bf3(r)
% bf3.m , p, mu are given parameters
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% define G(r)= F(r)-p = {integral of binormal pdf over [0,r]} - p
global p mu
y=quad('bf1',0,r) - p;
