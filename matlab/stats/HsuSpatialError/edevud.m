function p=edevud(abc,sxyz)
% edevud.m  2-21-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Ellipsoidal PDF, uncorrelated, ellipsoidal volume
%  abc=[ a, b, c]
% sxyz=[sx,sy,sz]

if min(abc)<=0 | min(sxyz)<=0
   error('check input arguments in edevud.m')
end

sxyz= sxyz./abc;
p=edsvud(1,sxyz);
