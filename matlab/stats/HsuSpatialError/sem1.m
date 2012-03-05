function p=sem1(r,s1,s2,ang)
% sem1.m  4-18-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% from Burt, single ellipse method #1
% ang = cut angle in degrees

d2r=pi/180 ;
% fictitious sigma and fictitious cut angle
  [sf, af, fac]=fsca(s1,s2,ang);

sa= sf/sqrt(2)/sin(af*d2r/2) ;
sb= sf/sqrt(2)/cos(af*d2r/2) ;
p=edcaud(r,sa,sb);
