%% Written by Andrew D. Wiles, 2009-12-03.
close all;
clear all;

%% parameters.
winsize = 100; % size of the window.
regmethod = 1;
toolNum = 1;
filename = '1000frames_12_4_001.csv';

%% TODO: get the rigid body definition.
refmrkcell = cell(1,2);
refmrkcell{1} = [
    -14.99      147.23       -0.01
    -0.00       69.65       -0.05
    0.00       -0.06        0.05
    -0.00        0.06      -60.07];
refmrkcell{2} = [
    0.00      134.92        0.00
    0.89      134.72      109.85
    0.00       -0.00        0.00
    0.00      -11.05       74.23];
refmrk = refmrkcell{toolNum};
tip = [0 0 0];
nMrks = size(refmrk,1);
[U, L, V0, x, A, Ainv] = xfrmToPA(refmrk);
for j =1:nMrks
    R0{j} = chol(A{j});
    [L0{j} U0{j}] = lu(A{j});
    
end

plotRigidBody(refmrk, tip);

%% initialize the data structures for each marker.
for j=1:nMrks
    fremean{j} = zeros(1,3);
    frecov{j} = zeros(3,3);
    frewin{j} = zeros(winsize, 3);
    sigmaFLE{j} = zeros(3,3);
end

%% TODO: get the parsed data.
data = readFromXLS(filename);
N = size(data{toolNum}{1},1);

frestack = zeros(nMrks*6,1);
Astack = [A{1};A{2};A{3};A{4}];
for i=1:N %N is the total Number of frames to parse.
    %% get the measured markers from the parsed data.
    for marks = 1:nMrks
        measmrk(marks,:) = data{toolNum}{marks}(i,:); % indexed by i.
    end
    [xfrm, newfre] = updateRegistration(refmrk, measmrk, regmethod);
    V = xfrm.rot*V0;
    for j=1:nMrks
        [fremean{j} frecov{j} frewin{j}] = updateFREStats(newfre(j,:), fremean{j}, frecov{j}, frewin{j});
        RPA{j} = V' * frecov{j} * V;
        %create FRE vec.
        frevec(1,1) = RPA{j}(1,1);
        frevec(2,1) = RPA{j}(1,2);
        frevec(3,1) = RPA{j}(1,3);
        frevec(4,1) = RPA{j}(2,2);
        frevec(5,1) = RPA{j}(2,3);
        frevec(6,1) = RPA{j}(3,3);

        frestack(6*(j-1)+1:6*j,1) = frevec;
        % estimate FLE using A.
        %flevec = Ainv{j}*frevec;
        %flevec = linsolve(A{j}, frevec);
        %flevec = A{j}\frevec;
        %flevec = R0{j} \ (R0{j}' \ frevec);
        flevec = inv(U0{j})*inv(L0{j})*frevec;
        %         if(j == 4)
        %             [flevec' frevec']
        %         end

        %build FLE covariance matrix
        fleest(1,1) = flevec(1,1);
        fleest(1,2) = flevec(2,1);
        fleest(2,1) = flevec(2,1);
        fleest(1,3) = flevec(3,1);
        fleest(3,1) = flevec(3,1);
        fleest(2,2) = flevec(4,1);
        fleest(2,3) = flevec(5,1);
        fleest(3,2) = flevec(5,1);
        fleest(3,3) = flevec(6,1);

        %rotate back into global axes.
        sigmaFLE{j} = V * fleest * V';
        flestats.rms2(i,j) = trace(sigmaFLE{j});
        flestats.rms(i,j) = sqrt(trace(sigmaFLE{j}));
        flestats.sigma11(i,j) = flevec(1,1);
        flestats.sigma12(i,j) = flevec(2,1);
        flestats.sigma13(i,j) = flevec(3,1);
        flestats.sigma22(i,j) = flevec(4,1);
        flestats.sigma23(i,j) = flevec(5,1);
        flestats.sigma33(i,j) = flevec(6,1);

        frestats.rms2(i,j) = trace(RPA{j});
        frestats.rms(i,j) = sqrt(trace(RPA{j}));
        frestats.sigma11(i,j) = frevec(1,1);
        frestats.sigma12(i,j) = frevec(2,1);
        frestats.sigma13(i,j) = frevec(3,1);
        frestats.sigma22(i,j) = frevec(4,1);
        frestats.sigma23(i,j) = frevec(5,1);
        frestats.sigma33(i,j) = frevec(6,1);
    end
    flestack = Astack\frestack;
    %build FLE covariance matrix
    fleest(1,1) = flestack(1,1);
    fleest(1,2) = flestack(2,1);
    fleest(2,1) = flestack(2,1);
    fleest(1,3) = flestack(3,1);
    fleest(3,1) = flestack(3,1);
    fleest(2,2) = flestack(4,1);
    fleest(2,3) = flestack(5,1);
    fleest(3,2) = flestack(5,1);
    fleest(3,3) = flestack(6,1);

    %rotate back into global axes.
    FLEStack = V * fleest * V';
    flestackstats.rms2(i,1) = trace(FLEStack);
    flestackstats.rms(i,1) = sqrt(trace(FLEStack));
    flestackstats.sigma11(i,1) = flestack(1,1);
    flestackstats.sigma12(i,1) = flestack(2,1);
    flestackstats.sigma13(i,1) = flestack(3,1);
    flestackstats.sigma22(i,1) = flestack(4,1);
    flestackstats.sigma23(i,1) = flestack(5,1);
    flestackstats.sigma33(i,1) = flestack(6,1);
    treIso(i,:) = (xfrm.rot * tip')' + xfrm.pos;
    %% Ramya's registration
    for j = 1:nMrks
        [V,D] = eig(sigmaFLE{j});
        W1 = V * inv(sqrt(D)) * V';
        W(:,:,j) = W1;
    end
    [R,t,fre,nIter,isConverged] = anisotropic_point_register(refmrk',measmrk',W,0.000001);
    xfrm.rot = R;
    xfrm.pos = t';
    newfreAniso = (xfrm.rot * refmrk')' + repmat(xfrm.pos, nMrks, 1) - measmrk;
    treAniso(i,:) = (xfrm.rot * tip')' + xfrm.pos;
    %% end of Ramya's registration
end

figure
hold on;
plot(1:N, flestats.rms(:,1));
plot(1:N, flestats.rms(:,2),'g');
plot(1:N, flestats.rms(:,3),'r');
plot(1:N, flestats.rms(:,4),'k');
plot(1:N, flestackstats.rms(:,1), 'm:');
hold off;
title('FLE RMS');
legend('Mrk 1', 'Mrk 2', 'Mrk 3', 'Mrk 4', 'Location', 'Best');

figure
subplot(3,2,1);
%plot(1:N, flestats.sigma11(:,1));
hold on;
%plot(1:N, flestats.sigma11(:,2),'g');
%plot(1:N, flestats.sigma11(:,3),'r');
plot(1:N, flestats.sigma11(:,4),'k');
%plot(1:N, flestackstats.sigma11(:,1), 'm:');
hold off;
title('FLE {\sigma_{11}}');
subplot(3,2,3);
%plot(1:N, flestats.sigma22(:,1));
hold on;
%plot(1:N, flestats.sigma22(:,2),'g');
%plot(1:N, flestats.sigma22(:,3),'r');
plot(1:N, flestats.sigma22(:,4),'k');
%plot(1:N, flestackstats.sigma22(:,1), 'm:');
hold off;
title('FLE {\sigma_{22}}');
subplot(3,2,5);
%plot(1:N, flestats.sigma33(:,1));
hold on;
%plot(1:N, flestats.sigma33(:,2),'g');
%plot(1:N, flestats.sigma33(:,3),'r');
plot(1:N, flestats.sigma33(:,4),'k');
%plot(1:N, flestackstats.sigma33(:,1), 'm:');
hold off;
title('FLE {\sigma_{33}}');
%legend('Mrk 1', 'Mrk 2', 'Mrk 3', 'Mrk 4', 'Overdetermined', 'Location', 'Best');

figure
hold on;
plot(1:N, frestats.rms(:,1));
plot(1:N, frestats.rms(:,2),'g');
plot(1:N, frestats.rms(:,3),'r');
plot(1:N, frestats.rms(:,4),'k');
hold off;
title('FRE RMS');
legend('Mrk 1', 'Mrk 2', 'Mrk 3', 'Mrk 4', 'Location', 'Best');

figure
subplot(3,2,1);
%plot(1:N, frestats.sigma11(:,1));
hold on;
%plot(1:N, frestats.sigma11(:,2),'g');
%plot(1:N, frestats.sigma11(:,3),'r');
plot(1:N, frestats.sigma11(:,4),'k');
%plot(1:N, frestackstats.sigma11(:,1), 'm:');
hold off;
title('FRE {\sigma_{11}}');
subplot(3,2,3);
%plot(1:N, frestats.sigma22(:,1));
hold on;
%plot(1:N, frestats.sigma22(:,2),'g');
%plot(1:N, frestats.sigma22(:,3),'r');
plot(1:N, frestats.sigma22(:,4),'k');
%plot(1:N, frestackstats.sigma22(:,1), 'm:');
hold off;
title('FRE {\sigma_{22}}');
subplot(3,2,5);
%plot(1:N, frestats.sigma33(:,1));
hold on;
%plot(1:N, frestats.sigma33(:,2),'g');
%plot(1:N, frestats.sigma33(:,3),'r');
plot(1:N, frestats.sigma33(:,4),'k');
%plot(1:N, frestackstats.sigma33(:,1), 'm:');
hold off;
title('FRE {\sigma_{33}}');

%figure
%plot(1:N, treIso(:,1)-treAniso(:,1));
%figure
%plot(1:N, treIso(:,2)-treAniso(:,2));
%figure
%plot(1:N, treIso(:,3)-treAniso(:,3));
%figure
%plot(1:N, sqrt((treIso(:,3)-treAniso(:,3)).^2+(treIso(:,2)-treAniso(:,2)).^2+(treIso(:,1)-treAniso(:,1)).^2));

