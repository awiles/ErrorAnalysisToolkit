function R=noval1(p)
% noval1.m   3-11-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% multiple ellipses with random axes orientations
% p=input( 'p=   (.5 in p.40 Burt) ');

d2r=pi/180;

load c:\mfile\nsd1.dat

nsd1
x=nsd1(:,1)        ;
y=nsd1(:,2)        ;
th=nsd1(:,3)*d2r   ;

a=cos(th).^2./x.^2 + sin(th).^2./y.^2
b=cos(th).^2./y.^2 + sin(th).^2./x.^2
c=cos(th).*sin(th).*(1./y.^2 - 1./x.^2)

rho=c./sqrt(a.*b)
w2=1./((1-rho.^2).*a)
z2=1./((1-rho.^2).*b)

sw2f=sum(w2)
sz2f=sum(z2)
rhof=sum(rho .*sqrt(w2) .*sqrt(z2))/(sqrt(sw2f)*sqrt(sz2f))

sx2f=0.5*(sw2f+sz2f+sqrt((sw2f+sz2f)^2-4*sw2f*sz2f*(1-rhof^2)))
sy2f=0.5*(sw2f+sz2f-sqrt((sw2f+sz2f)^2-4*sw2f*sz2f*(1-rhof^2)))

thf=atan2(2*rhof*sqrt(sw2f*sz2f),sw2f-sz2f)/d2r/2

sa=sqrt(sx2f)
sb=sqrt(sy2f)
% -------------------------------------------------

mu=min(sa,sb)/max(sa,sb) ;

K=bf4(p,mu)                 ;
R=K*max(sa,sb)              ;

disp(['R=CEP(' num2str(p*100) ')= ' num2str(R)])
disp(['R= ' num2str(K) ' times max{ sigma(x),sigma(y) }' ])
disp(['sigma(x)= ' num2str(max(sa,sb)) ])
disp(['sigma(y)= ' num2str(min(sa,sb)) ])
disp(['u= ' num2str(mu) ])
