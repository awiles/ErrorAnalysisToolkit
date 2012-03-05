function r2p3d(k)
% r2p3d.m  2-26-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
%
% for 3-dim normal distribution f(x,y,z) with  sigma(x)=sigma(y)=sigma(z)
% given k = multiple of SIGMA
% given k = multiple of MRE
% given k = multiple of RMSR
% given k = multiple of SEP
% r =k * SIGMA
% r =k * MSRE         (MRE= Mean Spherical Radial Error)
% r =k * RMSR         (RMSR= Root Mean Square Spherical Radial Error)
% r =k * SEP          (SEP= Spherical Error Probable)
% find p=F(r)
% where F(r) = integral of f(x) from 0 to r.
%     Maxwell pdf used for f(x) here.
%

disp( '  ')
disp( '  ')
ic=menu('Options:','k*SIGMA','k*MRE','k*RMSR','k*SEP','all the above');

if ic==1
r= k;
p=sf2(r) ;
disp(['For R=' ,num2str(k), ' * SIGMA'])
disp(['p=', num2str(p)])
end

if ic==2
r= k* 2*sqrt(2/pi);
p=sf2(r)     ;
disp(['For R=' ,num2str(k), ' * MRE'])
disp(['p=', num2str(p)])
end

if ic==3
r= k*sqrt(3);
p=sf2(r)     ;
disp(['For R=' ,num2str(k), ' * RMSR'])
disp(['p=', num2str(p)])
end

if ic==4
r=k*sf4(0.5);
p=sf2(r) ;
disp(['For R=' ,num2str(k), ' * SEP'])
disp(['p=', num2str(p)])
end

if ic==5
r= k;
p=sf2(r)      ;
disp(['For R=' ,num2str(k), ' * SIGMA'])
disp(['    p=', num2str(p)])

r= k* 2*sqrt(2/pi);
p=sf2(r)        ;
disp(['For R=' ,num2str(k), ' * MRE'])
disp(['    p=', num2str(p)])

r= k*sqrt(3);
p=sf2(r)     ;
disp(['For R=' ,num2str(k), ' * RMSR'])
disp(['    p=', num2str(p)])

r=k*sf4(0.5);
p=sf2(r)       ;
disp(['For R=' ,num2str(k), ' * SEP'])
disp(['    p=', num2str(p)])
end
