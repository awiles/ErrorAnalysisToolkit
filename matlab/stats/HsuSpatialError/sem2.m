function p=sem2(r,sa,sb,ang)
% sem2.m  4-18-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Example from Burt   method 2
% given r, s1, s2 and cut angle find probability
% p=P(r,s1,s2,cut_angle)

d2r=pi/180;
sa2=sa^2;
sb2=sb^2;
f=sin(ang*d2r)^2;

sx2=(sa2+sb2+sqrt((sa2+sb2)^2-4*f*sa2*sb2))/(2*f);
sy2=(sa2+sb2-sqrt((sa2+sb2)^2-4*f*sa2*sb2))/(2*f);

sx=sqrt(sx2)                                     ;
sy=sqrt(sy2)                                     ;
mu=min(sx,sy)/max(sx,sy) ;
k=r/max(sx,sy);
p=bf2(k,mu)   ;
clc
disp(['P(R)= ' num2str(p*100) '%, for R= ' num2str(r)])
disp(['R=r*sigma(x)= ' num2str(k) ' times sigma(x)' ])
disp(['sigma(x)= ' num2str(max(sx,sy)) ])
disp(['sigma(y)= ' num2str(min(sx,sy)) ])
disp(['u= ' num2str(mu) ])
