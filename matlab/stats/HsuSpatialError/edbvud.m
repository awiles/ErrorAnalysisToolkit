function p=edbvud(abc,sxyz)
% edbvud.m   2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% The probability of a
% 3D joint ellipsoidal normal pdf with (sx, sy, sz), rho=0
% over a box volume:  |x| < a,   |y| < b,  |z| < c.
%  abc=[ a, b, c]
% sxyz=[sx,sy,sz]
%

a=abc(1);
b=abc(2);
c=abc(3);

sx=sxyz(1);
sy=sxyz(2);
sz=sxyz(3);

if sx<=0 | sy<=0 | sz<=0 | a<0 | b<0 | c <0
   error('check signs of input argument')
end

p=nf2(a/sx)*nf2(b/sy)*nf2(c/sz);
