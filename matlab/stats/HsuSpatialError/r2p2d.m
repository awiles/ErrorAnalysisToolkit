function r2p2d(k)
% r2p2d.m   6-07-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% for 2-dim normal distribution with sigma(x)=sigma(y), rho=0.
% given k = multiple of SIGMA
% given k = multiple of MRE
% given k = multiple of DRMS(RMSR)
% given k = multiple of CEP
% r =k * SIGMA
% r =k * MRE         (MRE= Mean Radial Error)
% r =k * DRMS(RMSR)  (RMSR=Root Mean Square Radial Error)
%                     DRMS=Distance Root Mean Square
% r =k * CEP         (CEP= Circular Error Probable)
% find p=F(r)
% where F(x) = integral of f(x) from 0 to r.
% Rayleigh pdf used for f(x) here.
%

disp( '  ')
disp( '  ')
ic=menu('Options:','k*SIGMA','k*MRE','k*DRMS','k*CEP','all the above');

if ic==1
r= k;
p=cf2(r) ;
disp(['For r=' ,num2str(k), ' * SIGMA'])
disp(['p=', num2str(p)])
end

if ic==2
r= k* sqrt(pi/2);
p=cf2(r)    ;
disp(['For r=' ,num2str(k), ' * MRE'])
disp(['p=', num2str(p)])
end

if ic==3
r= k* sqrt(2);
p=cf2(r)    ;
disp(['For r=' ,num2str(k), ' * DRMS(RMSR)'])
disp(['p=', num2str(p)])
end

if ic==4
r=k*cf3(0.5);
p=cf2(r) ;
disp(['For r=' ,num2str(k), ' * CEP'])
disp(['p=', num2str(p)])
end

if ic==5
r= k;
p=cf2(r)     ;
disp(['For r=' ,num2str(k), ' * SIGMA'])
disp(['    p=', num2str(p)])

r= k* sqrt(pi/2);
p=cf2(r)       ;
disp(['For r=' ,num2str(k), ' * MRE'])
disp(['    p=', num2str(p)])

r= k* sqrt(2);
p=cf2(r)    ;
disp(['For r=' ,num2str(k), ' * DRMS(RMSR)'])
disp(['    p=', num2str(p)])

r=k*cf3(0.5);
p=cf2(r)      ;
disp(['For r=' ,num2str(k), ' * CEP'])
disp(['    p=', num2str(p)])
end
