function pf3(fname)
% pf3.m   01-25-96
% print figure to an eps file
%
%  Spatial Error Analysis ToolBoxc Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
fname=['e:\temp\hsu_test\' fname '.eps'];
% edit by ADW -- March 27, 2006 for Matlab 6.5
s1=['print -deps ' fname];
eval(s1)
