% chi1.m  1-15-94
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
clear
ic=2;   % skip step 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ic==1
% part I
%        use y=sqrt(x/n), x has chi-square distribution with deg n
clear global n p
global  n p

q=[0.5 .8 .9 .95];
l=[(0.1:0.025:0.2) (0.3:.1:.9) (1:9) (10:10:100)]';
NL = length(l);
NQ = length(q);

for jj=1:NQ
  p=q(jj);
  x =1;      %initial guess for iteration
% x =0.01;   %initial guess for fsolve
for ii=1:NL
  n=l(ii);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use fsolve from optimization toolbox %%
%%  xs=fsolve('x2disc2',x);
%%  xo= xs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% iteration %%%%%%%%%%%%

      p0=x2cdf(x);
   for j=1:2:9
      while (p0 > p)
      x=x-(0.1)^j;
      p0=x2cdf(x);
      end

      while (p0 < p)
      x=x+(0.1)^j;
      p0=x2cdf(x);
      end
   end
     xo=x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  r(ii,jj)=sqrt(xo/n)
  pause
end
end

DAT =[l  r];
save \mfile\chi.mat  l r

fname='c:\mfile\chi.dat'
fid=fopen(fname,'wt')
fprintf(fid,'%8.2f %8.4f %8.4f %8.4f %8.4f\n',DAT')
fclose(fid)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ic==2
% part II
load c:\mfile\chi.mat
NQ =4
adj= [ 0.7  0.5  0.15  -0.1];

semilogx(l,r)
xlabel('n')
ylabel(' R(p)/RMS ')
% title(' R(p)/RMS vs n,  p=0.5, 0.8, 0.9, 0.95')
grid

for kk=1:NQ
if kk==1   p=0.5
   elseif kk==2 p=0.8
   elseif kk==3 p=0.9
   elseif kk==4 p=0.95
end

   text(.5,r(5,kk)+adj(kk),['p=' num2str(100*p) '%'])
end
axis('square')
figure(1)
pause
pf3('f62')
end
