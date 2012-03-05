function [sf, af, fac]=fsca(sa,sb,ang)
% fsca.m   3-11-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% p. 1219 Bowditch, Fig. 6qi
% Find fictitious sigma, sigma factor, and fictitious cut angle
% ang, af in degrees

d2r=pi/180;
u=min(sa,sb)/max(sa,sb);

b=atan(u);
 sf=sin(2*b).*sqrt(sa^2+sb^2)/sqrt(2);
 af=asin(sin(2*b)*sin(ang*d2r))/d2r;
 fac = sf/max(sa,sb);

 disp([' fictitious sigma = ' num2str( sf )])
 disp([' fictitious cut angle = ' num2str(af) ' degrees'] )
 disp([' sigma factor = ' num2str(fac)] )
