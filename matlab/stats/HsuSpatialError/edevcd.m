function p=edevcd(abc,sxyz,rho)
% edevcd.m  12-26-98
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Ellipsoidal PDF, correlated, ellipsoidal volume
% abc =[a, b, c]
% sxyz=[sx, sy, sz]
% rho =[r12, r23, r13]

if min(abc)<=0 | min(sxyz)<=0 | max(abs(rho)) > 1
   error('check input arguments in edevcd.m')
end

if rho== [0 0 0]
   p=edevud(abc,sxyz);
else

a=abc(1);
b=abc(2);
c=abc(3);

sx=sxyz(1)/a;
sy=sxyz(2)/b;
sz=sxyz(3)/c;

r12=rho(1);
r23=rho(2);
r13=rho(3);

f1=1/(1+2*r12*r13*r23-r12^2-r13^2-r23^2);

   A= f1*(1-r23^2)/sx^2  ;
   B= f1*(1-r13^2)/sy^2  ;
   C= f1*(1-r12^2)/sz^2  ;

   q12=-(r12-r13*r23)*f1/(sx*sy);
   q13=-(r13-r23*r12)*f1/(sx*sz);
   q23=-(r23-r12*r13)*f1/(sy*sz);

   Q=[A q12 q13; q12 B q23; q13 q23 C];

   [V,D]=eig(Q);
   s1=1/sqrt(D(1,1));
   s2=1/sqrt(D(2,2));
   s3=1/sqrt(D(3,3));
   p=edsvud(1,[s1,s2,s3]);

end
