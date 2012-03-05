function p=noval(r)
% noval.m   3-11-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% multiple ellipses all axes parallel
% r=input ('r= (20 in p.32 Burt) ')            ;
load c:\mfile\nsd.dat

x=nsd(:,1)
y=nsd(:,2)

sa=sqrt(sum(x.*x))
sb=sqrt(sum(y.*y))

ang=90;
% ang=input ('cut angle=  (90 in p.32 Burt) ')  ;

disp(' ')
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
disp(['p(r)= ' num2str(p*100) '%, for r= ' num2str(r)])
disp(['r= ' num2str(k) ' times sigma(x)' ])
disp(['sigma(x)= ' num2str(max(sx,sy)) ])
disp(['sigma(y)= ' num2str(min(sx,sy)) ])
disp(['u= ' num2str(mu) ])
