%  tandn.m    3-23-96    plot Figure 2.3
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
clear global n
global n
n=100;
subplot(211)
fplot('nf1',[-5,5])
xlabel('x')
ylabel('f(x)')
title('N(0,1) PDF')
axis([-5 , 5, 0, 0.5])
grid

n=100;
subplot(212)
fplot('tdis',[-5,5])
grid
axis([-5 , 5, 0, 0.5])
sxlabel('\tau')
sylabel('g(\tau{})')
stitle('Student-t_n PDF, n=100')
fixstext
figure(1)
pf3('f23')
