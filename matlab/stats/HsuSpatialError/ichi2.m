function x=ichi2(p,N)
% ichi2.m   9-29-97
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
   s=['c:\mfile\p2x' num2str(N)];
   s1=['load ' s];
   eval(s1)
   x = interp1(pp,xx,p);
