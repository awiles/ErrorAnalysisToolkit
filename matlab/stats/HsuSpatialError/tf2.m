function p=tf2(r,u,v)
% tf2.m   2-20-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% for triple integration over [0,r], u, v as parameters
% avoid u=0 and v=0  in tf1.m

c0=-1/2;
c1=(v^2-1)/(2*v^2);
c2=(u^2-1)/(2*u^2);
fac1=4/(sqrt(2*pi)*(pi*u*v));

mr2=32;    % 32
mw2=64 ;   % 64
mt2=36;    % 36

dr=r/mr2;
dw=1/mw2;
dt=(pi/2)/mt2;

%           w.r.t.  r
for ir=1:mr2    % skip ir=0
r=dr*ir;

%           w.r.t.  w
       for iw=0:mw2
       w=dw*iw;

%        w.r.t.  t
          if c2==0 | w==1
           sum3=pi/2;
          else
           it=0:mt2;
           t=dt*it;
           y3=exp(c2*r.^2.*(1-w.^2).*cos(t).^2);
           sum3=sum(y3(1:mt2:mt2+1))+2*sum(y3(3:2:mt2-1))+4*sum(y3(2:2:mt2));
           sum3=sum3*dt/3;
          end

       if c1==0 | w==0
       y2(iw+1)=sum3;
       else
       y2(iw+1)=exp(c1*r.^2.*w.^2)*sum3  ;
       end

       end

sum2=sum(y2(1:mw2:mw2+1))+2*sum(y2(3:2:mw2-1))+4*sum(y2(2:2:mw2));
sum2=sum2*dw/3;
y1(ir+1)=exp(c0*r.^2).*r.^2*sum2;
end

sum1=y1(mr2+1)+2*sum(y1(3:2:mr2-1))+4*sum(y1(2:2:mr2));
p=fac1*sum1*dr/3;

if p>1
  p=1;
end
