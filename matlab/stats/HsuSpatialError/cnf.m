function p=cnf(ep,n)
% cnf.m  7-01-94
% determine confidence coeff about s, the sample estmate of population
% sigma for a given sample size and variation epsilon (in % of s)
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%

%            n       _
%      s^2 =SUM(xi - x)^2/(n-1)  is the sample estimate for sigma^2
%           i=1

%without the next line, n will stay the same after the first use.
clear global n

% n is used in the definition of chi-square CDF: x2cdf below.
global n

n=n-1;

if n<0 | ep<0
disp(' *** error in input *** ')
else

   if n<101 & n>=0                     %   small sample
       if n==0
         p=0;
       else
%                               s(1-epsilon) < sigma < s(1+epsilon)
         a=n/(1+ep)^2 ;
         b=n/(1-ep)^2 ;
         p=x2cdf(b)-x2cdf(a);
       end
   else                               %   large sample
         a=sqrt(2*n)/(1+ep)-sqrt(2*n-1);
         b=sqrt(2*n)/(1-ep)-sqrt(2*n-1);
         p=quad8('nf1',a,b);
   end

end
