function x=newton(xo,fun,dfun,tol)
% newton.m   1-1-97
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% xo is initial guess
% fun and dfun are string names

% example:
%         fun1 = sf3a;
%        dfun1 = sf1 ;
%      x = newton(xo,'sf3a','sf1',0.0001);

global p

  eval(['x = xo - ', fun,'(xo)/', dfun,'(xo);']);
  d=abs(x-xo);
  N=1;

  while d>tol
     xo=x ;
     eval(['x = xo - ', fun,'(xo)/', dfun,'(xo);']);
     N=N+1;
     d=abs(x-xo);
  end
