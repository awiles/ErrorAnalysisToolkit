function p2r2d(p)
% p2r2d.m  6-07-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given p and f(x), find r such that F(r)-p = 0,
% express r in terms of sigma, CEP, MRE, DRMS(RMSR)
% where F(r) = integral of f(x) from 0 to r.
%      Rayleigh pdf used for f(x) here.
%

disp('   ')
disp('   ')
disp([' For p=', num2str(p) ])
ic=menu('options:', 'k*sigma', 'k*MRE',...
                   'k*CEP', 'k*DRMS(RMSR)', 'all 4 above');
r0=cf3(p) ;
%

if ic==1
disp(['r=', num2str(r0), ' * SIGMA'])
end

if ic==2
r1=r0/sqrt(pi/2);
disp(['r=', num2str(r1), ' * MRE'])
end

if ic==3
r2=r0/cf3(.5);
disp(['r=', num2str(r2), ' * CEP'])
end

if ic==4
r3=r0/sqrt(2);
disp(['r=', num2str(r3), ' * DRMS(RMSR)'])
end

if ic==5
disp(['r=', num2str(r0), ' * SIGMA'])

r1=r0/sqrt(pi/2);
disp(['r=', num2str(r1), ' * MRE'])

r2=r0/sqrt(2);
disp(['r=', num2str(r2), ' * DRMS(RMSR)'])

r3=r0/cf3(0.5);
disp(['r=', num2str(r3), ' * CEP'])
end
