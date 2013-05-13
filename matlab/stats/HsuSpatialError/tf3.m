function r=tf3(p,u,v)
% tf3.m,  2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% solve for r from F(r,u,v) - p=0,
% sig(x)~=sig(y)~=sig(z)

if p>1 | p<0  | u<0 | v<0 | u>1 | v>1 | v>u
  error( 'input error for tf3(p)')
end

if p==1
    r=Inf;
else
    r=1;
    p0=tf2(r,u,v);
%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
 for j=1:2:5
     while (p0 > p)
       r=r-(0.1)^j;
       p0=tf2(r,u,v);
     end

     while (p0 < p)
       r=r+(0.1)^j;
       p0=tf2(r,u,v);
     end
 end
end
