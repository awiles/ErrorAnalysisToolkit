function R=edsvci(p,sxyz,rho)
% edsvci.m  2-27-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Ellipsoidal PDF, correlated, spherical volume
% sxyz=[sx,sy,sz]
% rho=[rxy,ryz,rxz]

if p<0 | p>1 | min(sxyz) < 0 | sxyz==[0 0 0] | min(rho)<-1 | max(rho)>1
   error('check input arguments in edsvci.m')
end

disp('Executing ... Patience')

R=mean(sxyz)*sf3(p);
p0 = edsvcd(R,sxyz,rho)  ;
%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:1:5
  while (p0 > p)
    R=R-(0.1)^j          ;
p0 = edsvcd(R,sxyz,rho)  ;
  end

  while (p0 < p)
    R=R+(0.1)^j          ;
p0 = edsvcd(R,sxyz,rho)  ;
  end
end
