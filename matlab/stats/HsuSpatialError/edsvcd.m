function p=edsvcd(R,sxyz,rho)
% edsvcd.m  12-26-98
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Ellipsoidal PDF, correlated, spherical volume
% sxyz=[sx, sy, sz]
% rho =[r12, r23, r13]

if R<0 | min(sxyz)<=0 | abs(max(rho)) > 1
   error('check input arguments in edsvcd.m')
end

sx=sxyz(1);
sy=sxyz(2);
sz=sxyz(3);

r12=rho(1);
r23=rho(2);
r13=rho(3);

f1=1/(1+2*r12*r13*r23-r12^2-r13^2-r23^2);

if rho==[0 0 0]
   p=edsvud(R,[sx,sy,sz]);
else
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
   p=edsvud(R,[s1,s2,s3]);
end
