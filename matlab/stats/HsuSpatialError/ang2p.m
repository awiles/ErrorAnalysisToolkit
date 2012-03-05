function p=ang2p(ang,p0)
% ang2p.m
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%

ro=sqrt(-2*log(1-p0))
d2r=pi/180;

    ang=ang*d2r;
    s1=csc(ang/2)/sqrt(2);
    s2=sec(ang/2)/sqrt(2);
    sx=max(s1,s2)
    sy=min(s1,s2);
    u=sy/sx
    mu=u;
    r=ro/sx
    p=bf2(r,mu)*100;
