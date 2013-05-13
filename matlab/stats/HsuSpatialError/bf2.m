function y=bf2(r,mu)
% bf2.m   1-15-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% Find probability  p
% function p=bf2(r,mu)
% where F(r) = integral of f(r,mu) over [0,r]
clear global mu
global mu

if r<0 | mu>1 | mu<0
        disp( ' input error for bf2(r,mu)' )
        disp( '   r >= 0 ' )
        disp( ' 0<= mu <= 1' )

elseif mu==0
        y=nf2(r);

elseif mu==1
        y=cf2(r);

else
        y=quad('bf1',0,r);
end
