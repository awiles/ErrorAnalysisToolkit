%  xandn.m    9-26-97    plot Figure 2.6
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
clear global n
global n
clf

subplot(211)
fplot('nf1',[-5,5])
sylabel('f(x)')
sxlabel('x')
stitle(['f(x)=N(0,1) is defined for all real u'])
axis([-5 , 5, 0, 0.5])
grid

n=100;
umin = -sqrt(2*n-1)
subplot(212)
fplot('gofu',[-5,5])
sylabel('g(u)')
sxlabel('u')
stitle(['g(u) is defined for u \geq ' num2str(umin), ', when n=100'])
axis([-5,5, 0, 0.5])
grid
figure(1)
 pf3('f26')
