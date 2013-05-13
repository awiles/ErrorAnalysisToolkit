function [gm , rms, ratio]=gmrms(fname)
% gmrms.m   6-25-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
xo=myfun1(fname)

[m,n]=size(xo);
if n==1
 x=xo ;
end

if n==2
 x=xo(:,1).^2 + xo(:,2).^2 ;
 x=sqrt(x);
end

if n==3
 x=xo(:,1).^2 + xo(:,2).^2  + xo(:,3).^2;
 x=sqrt(x);
end

 tmp=prod(x);
 gm=tmp^(1/m);

 x2=x.^2;
 tmp2=sum(x2);
 rms=sqrt(tmp2/m);

 ratio=gm/rms;
