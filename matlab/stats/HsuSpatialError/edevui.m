function [a,b,c]=edevui(p,k1,k2,sxyz)
%edevui.m  2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% ellipsoidal normal pdf, x,y,z are uncorrelated
% ellipsoidal volume of integration
%    x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
%
%  k1=a/b, k2=b/c
%
sx=sxyz(1);
sy=sxyz(2);
sz=sxyz(3);

if sx< 0 | sy< 0 | sz <0 | k1 <0 | k2 <0 |sxyz==[0 0 0]
 error('check input arguments in edevui.m')
end

disp('Executing ... Patience')

b=sf3(p)*mean(sxyz);
a=k1*b;
c=b/k2;
abc=[a b c]         ;
p0=edevud(abc,sxyz) ;

%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:1:5
  while (p0 > p)
    b=b-(0.1)^j;
    a=k1*b;
    c=b/k2 ;
    abc=[a b c]     ;
p0=edevud(abc,sxyz) ;
  end

  while (p0 < p)
    b=b+(0.1)^j;
    a=k1*b ;
    c=b/k2;
    abc=[a b c]     ;
p0=edevud(abc,sxyz) ;
  end
end
