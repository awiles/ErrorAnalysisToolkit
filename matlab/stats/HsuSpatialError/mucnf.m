function p=mucnf(ep,n)
% mucnf.m  3-07-94
% determine confidence level about xb, the sample estmate of population
% mean mu, for a given sample size and variation epislon (in % of s).
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%

%without the next line, n will stay the same after the first use.
clear global n

% n is used in the definition of student-t distribution below.
global n

capn=n;
n=n-1;

if n<0 | ep<0
disp(' *** error in input *** ')
else

   b=sqrt(capn)*ep ;

   if n<101 & n>=0                     %   small sample
       if n==0
         p=0;
       else
%                               xb-Äs  < å < xb+Äs     xb= x bar
         p=2*quad('tdis',0,b);
       end

   else                               %   large sample

       if b>4
         p=1;
       else
         p=2*quad('nf1',0,b);
       end
   end

end
