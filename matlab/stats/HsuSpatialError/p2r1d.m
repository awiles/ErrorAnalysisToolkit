function p2r1d(p)
% p2r1d.m  5-04-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% 1-dim normal distribution
% given p and f(x), find r such that F(r)-p = 0,
% express r in terms of SIGMA, MAE, LEP
% where F(r) = integral of f(x) from -r to r.
%   standard noraml distribution N(0,1) used for f(x) here.
%
r=nf3(p);

disp('   ')
disp('   ')
disp([' For p=', num2str(p) ])
ic=menu('options:', 'k*SIGMA', 'k*MAE', 'k*LEP', 'all the above');

if ic==1
disp(['r=', num2str(r), ' * SIGMA'])
end

if ic==2
r=r*sqrt(pi/2);
disp(['r=', num2str(r), ' * MAE'])
end

if ic==3
r=r/nf3(0.5);
disp(['r=', num2str(r), ' * LEP'])
end

if ic==4
disp(['r=', num2str(r), ' * SIGMA'])

r=r*sqrt(pi/2);
disp(['r=', num2str(r), ' * MAE'])

r=nf3(p)/nf3(0.5);
disp(['r=', num2str(r), ' * LEP'])
end
