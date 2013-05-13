function p=nf2a(r)
% nf2a.m   5-4-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given r and f(x), find p
% p = F(r) = integral of f(x) from -r to r
% f(x) is 1-dim normal distribution N(0,1)   with graph

clf

disp([' area under pdf between ' num2str(-r) ' and ' num2str(r)])

p=nf2(r);

b=abs(r) ; a=-abs(r);

fplot('nf1',[-4,4])
hold

m=10;
dx=(b-a)/m;

   for k=1:m-1
       x=a+k*dx;
       y=nf1(x);
       plot([x,x],[0,y],':')
   end
   ya=nf1(a);
   yb=nf1(b);
   plot([a,a],[0,ya],'--')
   plot([b,b],[0,yb],'--')

xlabel(' x')
ylabel(' f(x) ')
t=['Area under Normal PDF f(x), for x in [' num2str(-r) ',',...
   num2str(r) '] is ' num2str(p,4) ];
title(t)
figure(1)
pf3('f11')   % or pf3('f33')
