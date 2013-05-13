function y=psin(x)
% psin.m  3-14-97
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% evaluate psi(x)
% psin(x)=psinh(2x)
% for x=integers, integer+0.5,  and 0<x<1
% recursive formula used for other values of x

eu=.5772156649;
global n
n=x;

 if n<0
   y=psin(1-x)-pi/tan(pi*x);
 else
 if n<1                               % n is less than 1
   y=quad8('psi',0,1-eps)-eu-1/n;     % 1-eps to avoid singularity
%  n: integer                         % n is integer
 elseif n-floor(n)==0 & n>=1
     k=1:n-1 ;  s=1./k;
     y= -eu + sum(s) ;
%  n: integer+0.5                     % n is integer+0.5
 elseif n-floor(n)==1/2  & n>1
     m=(2*n-1)/2;
     k=1:m ;  s=1./(2*k-1);
     y= -eu -2*log(2) + 2*sum(s) ;
 else
     x=n-1;
     y=psin(x)+ 1/x;
 end
 end

clear global n
