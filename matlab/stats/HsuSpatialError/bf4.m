function r=bf4(p,mu)
% bf4.m  2-14-94
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% find r for G(r)=F(r)-p=0

%clear global p mu
%global p mu

if p>1 | p<0 | mu>1 | mu < 0
   disp( ' input error for bf4(p,mu) ')
   disp( ' 0 <= p  <= 1 ')
   disp( ' 0 <= mu <= 1 ')

elseif mu==1
   r=cf3(p);

elseif mu==0
   r=nf3(p);

else
   if p==1
     r=Inf;
   else
t1=clock;
%    r=fsolve('bf3',1);
     r=newton(1,'bf3','bf1',0.0001);
time_used=etime(clock,t1)
   end
end
