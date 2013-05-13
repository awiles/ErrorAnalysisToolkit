function R=edsvui(p,sxyz)
% edsvui.m  2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Ellipsoidal PDF, uncorrelated, spherical volume
% sxyz=[sx,sy,sz]

disp(' Executing ...  Patience ')
if p<0 | p>1 | min(sxyz) < 0 | sxyz==[0 0 0]
   error('check input arguments in edsvui.m')
end

R=mean(sxyz)*sf3(p);
p0 = edsvud(R,sxyz);
%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:1:5
  while (p0 > p)
    R=R-(0.1)^j    ;
p0 = edsvud(R,sxyz);
  end

  while (p0 < p)
    R=R+(0.1)^j    ;
p0 = edsvud(R,sxyz);
  end
end
