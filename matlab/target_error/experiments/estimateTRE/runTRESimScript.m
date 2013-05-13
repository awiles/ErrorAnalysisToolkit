% common parameters.
parm.experimentDir = 'E:\docs\research\phd\experiments\TRE Experiments';
parm.M = 1000;
parm.N =10000;
fleSigmaType = 'aniso';
fleSigma = diag([0.0995 0.0995 0.2985].^2);
toolDesign = 'e';
A = 71;
B = 0.5*54;
rho = 85*[1:8];


for i = 1:8
    parm.case{i}.fleSigma = fleSigma;
    [parm.case{i}.Rigid.mrk, parm.case{i}.Rigid.normals, parm.case{i}.Rigid.tip]...
        = getWestToolDesign(toolDesign, A, B, rho(i));
    parm.case{i}.xfrm.pos = [0 0 -1000];
    parm.case{i}.xfrm.R = eye(3);
    parm.case{i}.varyRotation = 1;
    parm.case{i}.varyRotationAngle = 30;
    parm.case{i}.key = rho(i);
    parm.case{i}.name = sprintf('Des_%s_%s_A%dB%drho%d',toolDesign, fleSigmaType,A,B,rho(i));
end

test  = runTRESimulations(parm);

% common parameters.
parm.experimentDir = 'E:\docs\research\phd\experiments\TRE Experiments';
parm.M = 1000;
parm.N =10000;
fleSigmaType = 'aniso';
fleSigma = diag([0.0995 0.0995 0.2985].^2);
toolDesign = 'e';
A = 0.5*71;
B = 0.5*54;
rho = 85*[1:8];


for i = 1:8
    parm.case{i}.fleSigma = fleSigma;
    [parm.case{i}.Rigid.mrk, parm.case{i}.Rigid.normals, parm.case{i}.Rigid.tip]...
        = getWestToolDesign(toolDesign, A, B, rho(i));
    parm.case{i}.xfrm.pos = [0 0 -1000];
    parm.case{i}.xfrm.R = eye(3);
    parm.case{i}.varyRotation = 1;
    parm.case{i}.varyRotationAngle = 30;
    parm.case{i}.key = rho(i);
    parm.case{i}.name = sprintf('Des_%s_%s_A%dB%drho%d',toolDesign, fleSigmaType,A,B,rho(i));
end

test  = runTRESimulations(parm);