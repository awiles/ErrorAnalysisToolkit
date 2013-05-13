% chi2.m  1-16-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% R(p)/RMS=r  vs  GM/RMS=z

clear; clf;

% result from psi1.m
load c:\mfile\psi.dat
z=psi(:,2);

% result from chi1.m
load c:\mfile\chi.dat
r=chi(:,2:5);

DAT=[z r]
save c:\mfile\psichi.mat z r


fname='c:\mfile\psichi.dat'
fid=fopen(fname,'wt')
fprintf(fid,'%8.4f%8.4f%8.4f%8.4f%8.4f\n',DAT')
fclose(fid)

plot(z,r)
grid

adj = [0 0 0.06 0.02];
for kk=1:4
if kk==1   p=0.5      ;
   elseif kk==2 p=0.8 ;
   elseif kk==3 p=0.9 ;
   elseif kk==4 p=0.95;
end

 text(0.2,r(7,kk)+adj(kk),['p=' num2str(100*p) '%'])
end

xlabel(' x=GM/RMS  ')
ylabel(' y=R(p)/RMS  ')
% title(' R(p)/RMS vs z=GM/RMS,  for p=0.5, 0.8, 0.9, 0.95')
axis('square')
figure(1)
pause
pf3('f63')
