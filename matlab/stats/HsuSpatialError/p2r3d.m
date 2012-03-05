function  p2r3d(p)
% p2r3d.m   2-26-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% r=r3ng(p)
% given p and f(x), find r such that F(r)-p = 0,
% express r in terms of SIGMA, MSRE, SEP, RMSR
% where F(r) = integral of f(x) from 0 to r.
%      Maxwell pdf used for f(x) here.
%
%

disp('   ')
disp('   ')
disp([' For p=', num2str(p) ])
ic=menu('options:', 'k*SIGMA', 'k*MRE', 'k*RMSR','k*SEP','all 4 above');
r=sf4(p);

if ic==1
disp(['R=', num2str(r), ' * SIGMA'])
end

if ic==2
r=r/sqrt(8/pi);
disp(['R=', num2str(r), ' * MRE'])
end

if ic==3
r=r/sqrt(3);
disp(['R=', num2str(r), ' * RMSR'])
end

if ic==4
r=r/sf4(0.5);
disp(['R=', num2str(r), ' * SEP'])
end

if ic==5
disp(['R=', num2str(r), ' * SIGMA'])

r1=r/sqrt(8/pi);
disp(['R=', num2str(r1), ' * MSRE'])

r2=r/sqrt(3);
disp(['R=', num2str(r2), ' * RMSR'])

r3=r/sf4(0.5);
disp(['R=', num2str(r3), ' * SEP'])
end
