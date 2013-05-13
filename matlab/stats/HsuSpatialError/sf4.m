function r=sf4(p)
% sf4.m  2-19-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% solve for r in G(r)=0
% where  G(r)= F(r)-p = {integral of trinormal pdf over [0,r]} - p
clear global p
global p
if  p>1  | p<0
   disp( ' error in input ')
   disp( ' 0  <= p  <= 1 ')
else
   if p==1
     r=Inf;
   elseif p==0
     r=0;
%  elseif p==0.025 | p==0.05
   elseif p<=0.05  & p>0
%    r=fsolve('sf3a',0.5);
     r=newton(0.5,'sf3a','sf1',0.0001);
   else
%    r=fsolve('sf3a',1);
     r=newton(1,'sf3a','sf1',0.0001);
   end
end
