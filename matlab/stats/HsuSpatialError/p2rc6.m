function R=p2rc6(gm,rms,p)
% p2rc6.m  1-2-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% when n= 0.9358
% gm=5.8054,  rms=11.5296,  p=0.5  ---> r=7.542
clear global n
global n

fid=fopen('c:\mfile\psi.dat','r');
dat=fscanf(fid,'%g %g',[2,inf]);
fclose(fid);

ratio=gm/rms;
n2=dat(1,:);
r2=dat(2,:);
n =interp1(r2,n2,ratio,'spline');

%%%%%% use fsolve.m (from optimization toolbox) %%%%%%%
% global p
% x=fsolve('x2disc2',1) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% iteration %%%%%%%%%%%%
if p==1
   x=Inf;
else
   x=1;  % initial guess
      p0= x2cdf(x);
   for j=1:2:5
      while (p0 > p)
      x=x-(0.1)^j;
      p0= x2cdf(x);
      end

      while (p0 < p)
      x=x+(0.1)^j;
      p0= x2cdf(x);
      end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   disp([' n = ', num2str(n)])
   R=rms*sqrt(x/n);
