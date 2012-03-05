%nf7.m
% nf7.m  3-06-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
function nf7(p)
clf

  load c:\mfile\eg3p2.mat

mu=mean(z)
sig=std(z)

   xx1=floor(mu-4*sig);
   xx2=ceil(mu+4*sig);
   d=(xx2-xx1)/100;
   x=[xx1:d:xx2];
   y=nf6(x,mu,sig);
   plot(x,y)
   hold

r=nf3(p)
a=sig*(-r)+mu;
b=sig*(r)+mu;

m=10;
dx=(b-a)/m;

   for k=1:m-1
       x=a+k*dx;
       y=nf6(x,mu,sig);
       plot([x,x],[0,y],':')
   end
   ya=nf6(a,mu,sig);
   yb=nf6(b,mu,sig);
   plot([a,a],[0,ya],'--')
   plot([b,b],[0,yb],'--')


xlabel(' z ')
ylabel(' g(z) ')
t=['Area under Normal PDF g(z), for z in [' num2str(a) ',' num2str(b) '] is ' num2str(p) ];
title(t)
figure(1)
pf3('f34')
