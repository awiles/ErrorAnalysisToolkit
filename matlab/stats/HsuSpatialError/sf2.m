function y=sf2(r)
% sf2.m  3-06-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% find p=F(r)
% where   F(r)={integral of Maxwell pdf from 0 to r}, in terms of erf(x)

if r<0
   disp( 'input error ')
   disp( '  r>=0      ')
else
fac1=sqrt(2/pi);
y=fac1*(-r.*exp(-r.^2/2))+erf(r/sqrt(2));
end
