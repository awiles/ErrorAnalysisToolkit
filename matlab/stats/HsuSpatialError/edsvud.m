function p=edsvud(R,sxyz)
% edsvud.m  2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Ellipsoidal PDF, uncorrelated, spherical volume
% sxyz=[sx,sy,sz]

if R<0 | min(sxyz) < 0 | sxyz==[0 0 0]
   error('check input arguments in edsvud.m')
end

sx=sxyz(1);
sy=sxyz(2);
sz=sxyz(3);

sig=sort(sxyz);
r=     R/sig(3);
u=sig(2)/sig(3);
v=sig(1)/sig(3);

p=edsv(r,u,v);
