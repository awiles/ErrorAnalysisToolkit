function p=r2pc6(gm,rms,r)
% r2pc6.m 6-22-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% when n= 0.9358  and r=14.7497
% gm=5.81, rms=11.53            ---> p=0.8000.
clear global n
global n

fid=fopen('c:\mfile\psi.dat','r');
dat=fscanf(fid,'%g %g',[2,inf]);
fclose(fid);

ratio=gm/rms;

n2=dat(1,:);
r2=dat(2,:);

n =interp1(r2,n2,ratio,'spline');
ms=rms^2;
x=n*r^2/ms;

    p=x2cdf(x) ;
