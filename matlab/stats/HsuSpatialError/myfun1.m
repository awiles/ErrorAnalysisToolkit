function x=myfun1(fname)
% myfun1.m  2-14-96
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
   s1=['load ' 'c:\mfile\' fname]
   eval(s1)
   eval(['x = ', fname(1:3) ';'])
