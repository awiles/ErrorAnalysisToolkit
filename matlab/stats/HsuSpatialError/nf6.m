function y=nf6(x,mu,sig)
% nf6.m  3-06-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% nf6.m  define f(x)=N(mean, variance)
% notice .^  is used  for vector x
fac1=1/(sqrt(2*pi)*sig);
y=fac1*exp(-(x-mu).^2/(2*sig^2));
