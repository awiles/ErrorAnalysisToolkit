clear all;

% define West tool design parameters.
A = 71;
B = 54;
rho = 170;

testnum = 0;
designs = {'d', 'e'};

for i = 1:length(designs)
    testnum = testnum +1;
    testname = sprintf('calcTRE.m Isotropic FLE, West Tool Design %s', designs{i});
    fprintf( 'Test #%d: %s started...\n', testnum, testname);
    % define the rigid body defintion in local probe coordinates.
    % West design in Fig. 2(e).
    [mrk, normals, tip] = getWestToolDesign(designs{i}, A, B, rho);
    
    % now that the probe is created, let's rotate the probe into a useful
    % orientation.
    rotZ = 0;
    R = getRotMatrixd([0, 0, rotZ]);
    mrk = (R * mrk')';
    tip = (R * tip')';
    
    % define an arbitrary FLE error value where each axis is identical.
    FLE = 0.33;
    
    % using West's formula, compute the TRE.
    [TRE_West, Sigma_West, SigmaPA_West] = calcTRE(FLE,[mrk;tip],'FigNum', 1, 'Isotropic');
    % using my anisotropic formula with isotropic --> is it the same?
    [TRE_Wiles, Sigma_Wiles, SigmaPA_Wiles] = calcTRE(FLE,[mrk;tip], 'FigNum', 2);
    
    if( abs(TRE_Wiles - TRE_West) < 0.0001)
        fprintf('\tTest #%d: %s passed.\n', testnum, testname);
    else
        fprintf('\tTest #%d: %s FAILED.\n', testnum, testname);
    end
    fprintf('\t\tTRE West  = %3.4f\n', TRE_West);
    fprintf('\t\tTRE_Wiles = %3.4f\n', TRE_Wiles);
end

%% let's look at the TRE with the RB in different orientations.
testnum = testnum +1;
testname = 'calcTRE.m Various FLE, West Tool Design e, Different Orientations';
fprintf( 'Test #%d: %s started...\n', testnum, testname);
ang = -45:1:45;
TRE_iso_iso = zeros(size(ang));
TRE_iso_ani = zeros(size(ang));
TRE_ani_iso = zeros(size(ang));
TRE_ani_ani = zeros(size(ang));

% define isotropic FLE.
FLE = 0.33;

for i = 1:length(ang)
    R = getRotMatrixd([0, ang(i), 0]);
    mrk_test = (R * mrk')';
    tip_test = (R * tip')';
    
    TRE_iso_iso(i) = calcTRE(FLE,[mrk_test;tip_test], 'Isotropic');
    TRE_iso_ani(i) = calcTRE(FLE, [mrk_test;tip_test] );
end

% define anisotropic FLE where z std. dev. is 3 times x and y.
FLE = [0.0995 0.0995 0.2985].^2;
%TRE_Wiles = calcTRE(FLE,[x;tip0]);

ang = [-45:1:45];

for i = 1:length(ang)
    R = getRotMatrixd([0, ang(i), 0]);
    mrk_test = (R * mrk')';
    tip_test = (R * tip')';
    
    TRE_ani_iso(i) = calcTRE(FLE,[mrk_test;tip_test], 'Isotropic');
    TRE_ani_ani(i) = calcTRE(FLE,[mrk_test;tip_test]);
end

figure(3);
plot(ang, TRE_iso_iso, 'b-.');
hold on;
plot(ang, TRE_iso_ani, 'g');
plot(ang, TRE_ani_iso, 'r--');
plot(ang, TRE_ani_ani, 'k:');
hold off;
xlabel('Rotation about y-axis from default orientation (degrees)');
ylabel('TRE (mm)');
legend('Isotropic \Sigma_{fle}, Isotropic Formula',...
    'Isotropic \Sigma_{fle}, Anisotropic Formula',...
    'Anisotropic \Sigma_{fle}, Isotropic Formula',...
    'Anisotropic \Sigma_{fle}, Anisotropic Formula');
title('Comparison of TRE Model using Isotropic and Anisotropic FLE and Formulas');
fprintf('\tTest #%d: %s passed.\n', testnum, testname);