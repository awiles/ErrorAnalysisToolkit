function y=pdft(td,ux,uy,sx,sy)
%  pdft.m  11-03-95
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
%  td in degrees
   sx2=sx^2; sy2=sy^2;
   d2r=pi/180;
   t=td*d2r;
   ct=cos(t);
   st=sin(t);
   ct2=ct^2;
   st2=st^2;
   A=ct2/sx2+st2/sy2;
   B=-2*(ux*ct/sx2+uy*st/sy2);
   C=(ux/sx)^2+(uy/sy)^2;
   w=B/sqrt(8*A);
   w2=w^2;
   y=1/(2*pi*sx*sy*A)*exp(-C/2)*(1-w*sqrt(pi)*erfc(w)*exp(w2));
