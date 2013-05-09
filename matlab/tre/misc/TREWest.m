% TREWest.m

% polaris active case 3 on pg. 541(9) of West.
A = 70;
B = 40;
rho = 170;
r = 63.6;
d = 400;

% probe marker and tip positions
xp = [-A/2 -B/2 0; -A/2 B/2 0; A/2 B/2 0; A/2 -B/2 0];
tp = [rho 0 0];

% coordinate reference frame
ang = deg2rad(45);
xc = [r*cos(ang) r*sin(ang) 0; -r*cos(ang) r*sin(ang) 0;...
        -r*cos(ang) -r*sin(ang) 0; r*cos(ang) -r*sin(ang) 0];
tc = [d 0 0];

% create a rotation matrix.
R = getRotMatrix(0, 55, 0);

xp = (R*xp')';
tp = (R*tp')';

FLE = 0.33;
%FLE = [0.0099 0.0099 0.0891];  %variances on the diagonal.
%stdev = sqrt(FLE^2/3);

TREp = calcTRE_NDI(FLE, [xp;tp],1);
fprintf('Probe TRE = %4.2f\n', TREp);

TREc = calcTRE_NDI(FLE, [xc;tc],2);
fprintf('CRF TRE = %4.2f\n', TREc);

TREt = sqrt(TREp^2 + TREc^2);
fprintf('Combined TRE = %4.2f\n', TREt);

