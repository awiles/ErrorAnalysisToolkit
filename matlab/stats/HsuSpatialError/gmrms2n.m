function n=gmrms2n(gm,rms)
% gmrms2n.m   3-14-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
fid=fopen('c:\mfile\psi.dat','r');
dat=fscanf(fid,'%g %g',[2,inf]);
fclose(fid);
ratio=gm/rms ;

n2=dat(1,:);
r2=dat(2,:);
n =interp1(r2,n2,ratio,'spline');
