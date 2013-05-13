function r2p1d(k)
% r2p1d.m  5-04-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given k = multiple of SIGMA
% given k = multiple of MAE
% given k = multiple of LEP
% r =k * SIGMA
% r =k * MAE         (MAE= Mean Absolute Error)
% r =k * LEP         (LEP= Linear Error Probable)
% find p = F(r) = integral of f(x) from -r to r
% f(x) is 1-dim normal distribution N(0,1)

disp( '  ')
disp( '  ')
ic=menu('Options:','k*SIGMA','k*MAE','k*LEP','all the above');

if ic==1
r= k;
p=nf2(r) ;
disp(['For r=' ,num2str(k), ' * SIGMA'])
disp(['p=', num2str(p)])
end

if ic==2
r= k* sqrt(2/pi);
p=nf2(r)     ;
disp(['For r=' ,num2str(k), ' * MAE'])
disp(['p=', num2str(p)])
end

if ic==3
r0=nf3(0.5);
r=k*r0;
p=nf2(r) ;
disp(['For r=' ,num2str(k), ' * LEP'])
disp(['p=', num2str(p)])
end

if ic==4
r= k;
p=nf2(r)      ;
disp(['For r=' ,num2str(k), ' * SIGMA'])
disp(['    p=', num2str(p)])

r= k* sqrt(2/pi);
p=nf2(r)        ;
disp(['For r=' ,num2str(k), ' * MAE'])
disp(['    p=', num2str(p)])

r=k*nf3(0.5);
p=nf2(r)       ;
disp(['For r=' ,num2str(k), ' * LEP'])
disp(['    p=', num2str(p)])
end
