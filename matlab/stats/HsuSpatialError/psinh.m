function y=psinh(n2)
% psinh.m   3-14-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% evaluate psi(n/2)
% psinh(2x)=psin(x)
% for n=integers and 0<n<1

eu=.5772156649;
 global n
 n=n2/2;
 m=floor(n);

% evaluate psi
 if n2<1 & n2 >0                 % n2 is less than 1
   y=quad8('psi',0,1-eps)-eu-1/n;       % 1-eps to avoid singularity
 end
% odd n2: 1,3,5,...,99
   if n2==m*2+1                  % n2 is odd
     k=1:m ;  s=1./(2*k-1);
     y=-eu - 2*log(2) + 2*sum(s) ;
   end
% even n2: 2,4,6,...,100
   if n2==m*2                    % n2 is odd
     k=1:m-1;  s=1./k;           % n2 is even
     y=-eu + sum(s) ;
   end
clear global n
% z=GM/RMS
% z=exp(y/2)/sqrt(n);
