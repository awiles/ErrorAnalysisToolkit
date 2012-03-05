% psi1.m   3-12-97
%
%  Spatial Error Analysis ToolBox Version 1.0,  October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% GM/RMS vs N
% part I
clf
eu=.5772156649;

x=[(.1:0.025:.2) (0.3:.1:.9) (1:9) (10:10:100)]   ;
Nx = length(x);

for j=1:Nx
 n2=x(j);    m=n2/2;
  y=psinh(n2);
% y=psin(m);
 z(j)=exp(y/2)/sqrt(m);
end

DAT=[x' z']
save c:\mfile\psi.mat x z

fname='c:\mfile\psi.dat'
fid=fopen(fname,'wt')
fprintf(fid,'%8.4f  %8.4f\n',DAT')
fclose(fid)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% part II
load c:\mfile\psi.mat

semilogx(x,z, x,z,'o')
axis('square')
xlabel('n')
ylabel('GM/RMS')
% title(' GM/RMS  vs. n ')
grid
figure(1)
pf3('f61')
