function estimateFLEStatic(casestudy, winsize)
%% simulator for estimating FLE.
N = 12000;
fle1 = diag([0.01, 0.01, 0.09]);
fle2 = 0.11*eye(3);
fle3 = diag([0.03, 0.01, 0.09]);
%winsize = 200;
covmethod = 1;
pamethod = 1; %method used to rotate into principal axes: 1: xfrmed positions, 2: true positions.
%% get the rigid body design
[refmrk, normals, tip] = getWestToolDesign('e', 71, 54, 85);
%refmrk = [refmrk; 100 20 0];
nMrks = size(refmrk,1);

[U, L, V0, x, A, Ainv] = xfrmToPA(refmrk);
A{1}
inv(A{1})
V0
%% orientation case studies.
%casestudy = 1;
%setxfrm.R = getRandOrientation();
switch casestudy
    case 1
        name = 'case1'
        setxfrm.R = eye(3);
        setxfrm.pos = [0 0 0];
    case 2
        name = 'case2'
        setxfrm.R = getRotMatrixd([5,0,30]);
        setxfrm.pos = [0 0 0];
end
%% run the Monte Carlo simulation.
mrk = (setxfrm.R * refmrk')' + repmat(setxfrm.pos, size(refmrk,1), 1);
%[RMS, SigmaCov, SigmaCovPA, A] = calcTRE(fle1, [mrk; tip] );
for i = 1:nMrks
    % Ainv{i} = inv(A{i});
    % perform some preallocations.
    frewin{i} = zeros(winsize,3);
    frecovwin{i} = zeros(winsize,6);
    fleValues{i} = zeros(6,N);
    fleEstValues{i} = zeros(6,N);
    fleEigValues{i} = zeros(3,N);
    fleEstEigValues{i} = zeros(3,N);
    fremean{i} = [0 0 0];
    frecov{i} = zeros(3);
    frecovmean{i} = [0 0 0 0 0 0]; %for the moving average.
    frenoisycov{i} = zeros(N,6);
    fresmoothcov{i} = zeros(N,6);
end
measfreRMS = zeros(N, nMrks);
fleRMS = zeros(N,1);
fleestRMS = zeros(N,nMrks);
meanFLERMS = zeros(N,1);
meanFLEEigValues = zeros(3,N);
%RRMS = zeros(N,nMrks);
frecovRMS = zeros(N,nMrks);
fresmoothRMS = zeros(N,nMrks);

%% set up the Kalman Filter parameters.
prevfrecov = [0 0 0 0 0 0]';
prevP = eye(6);
Q = eye(6);
R = Q;
%R = diag([0.1, 0.01, 0.01, 0.1, 0.01, 0.1]);
K = prevP*inv(prevP + R);

%% run the loop.
starttime = cputime;
for i = 1:N
    %set up the fle.
    if(i<=N/3)
        fle = fle1;
    elseif(i>N/3 && i < 2*N/3)
        fle = fle2;
    else
        fle = fle3;
    end
    %R = fle;
    %add the noise.
    measmrk = mrk + getMeasNoise(fle, nMrks);

    %*****************************************
    %compute the transform and xfrm positions.
    %*****************************************
    % using SVD:
    xfrm = getRigidXfrmSVD( refmrk, measmrk );
    xfrmmrk = (xfrm.rot * refmrk')' + repmat(xfrm.pos, nMrks, 1);
    
    % using Horn:
    %xfrm = getRigidXfrm( refmrk, measmrk );
    %xfrmmrk = (quat2rm(xfrm.rot) * refmrk')' + repmat(xfrm.pos, nMrks, 1);
    
    switch pamethod
        case 1
            %[U, L, V, x, A, Ainv] = xfrmToPA(xfrmmrk);
            %[U, L, V] = xfrmToPA(xfrmmrk);
            %[V, D] = getPA(xfrmmrk);
            V = xfrm.rot*V0;
            %V = V0*xfrm.rot';
        case 2
            %[U, L, V, x, A, Ainv] = xfrmToPA(mrk);
            %[U, L, V] = xfrmToPA(mrk);
            [V, D] = getPA(mrk);
        otherwise
            error('Invalid method to rotate into PA given.');
    end
    %compute the FRE.
    fre = xfrmmrk - measmrk;
    meanFLE = zeros(3);
    for j=1:nMrks
        switch covmethod
            case 1
                %update the mean.
                fremean{j} = fremean{j} + 1/winsize* (fre(j,:) - frewin{j}(end,:));
                %update the covariance
                frecov{j} = frecov{j} + 1/(winsize-1) * (fre(j,:)'*fre(j,:) - frewin{j}(end,:)'*frewin{j}(end,:));
                %update the window.
                frewin{j} = circshift(frewin{j},1);
                frewin{j}(1,:) = fre(j,:);
            case 2
                %update the window.
                frewin{j} = circshift(frewin{j},1);
                frewin{j}(1,:) = fre(j,:);
                fremean{j} = mean(frewin{j});
                frecov{j} = cov(frewin{j}); %TODO: replace this with more efficient update code.
            otherwise
                error('Invalid covariance method given.');
        end

        %rotate R into PA.
        %RPA{j} = V' * R{j} * V;
        RPA{j} = V' * frecov{j} * V;

        %create FRE vec.
        frevec(1,1) = RPA{j}(1,1);
        frevec(2,1) = RPA{j}(1,2);
        frevec(3,1) = RPA{j}(1,3);
        frevec(4,1) = RPA{j}(2,2);
        frevec(5,1) = RPA{j}(2,3);
        frevec(6,1) = RPA{j}(3,3);

        %try moving average.
        %update the mean.
        frecovmean{j} = frecovmean{j} + 1/winsize* (frevec' - frecovwin{j}(end,:));
        %update the window.
        frecovwin{j} = circshift(frecovwin{j},1);
        frecovwin{j}(1,:) = frevec';
        % store for plotting
        frenoisycov{j}(i,:) = frevec;
        fresmoothcov{j}(i,:) = frecovmean{j}';
        fresmoothRMS(i,j) = sqrt(frecovmean{j}(1) + frecovmean{j}(4) + frecovmean{j}(6));
        %estimate the FLE
        %flevec = Ainv{j}*frevec;
        flevec = Ainv{j}*frecovmean{j}';

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
        meanFLE = meanFLE + fleest;

        %TODO:rotate back into global axes.
        sigmaFLE = V * fleest * V';

        measfreRMS(i,j) = sqrt(fre(j,:)*fre(j,:)');
        %freRMS(i) = sqrt(curfre'*curfre);
        fleValues{j}(1,i) = fle(1,1);
        fleValues{j}(2,i) = fle(1,2);
        fleValues{j}(3,i) = fle(1,3);
        fleValues{j}(4,i) = fle(2,2);
        fleValues{j}(5,i) = fle(2,3);
        fleValues{j}(6,i) = fle(3,3);

        fleEstValues{j}(1,i) = sigmaFLE(1,1);
        fleEstValues{j}(2,i) = sigmaFLE(1,2);
        fleEstValues{j}(3,i) = sigmaFLE(1,3);
        fleEstValues{j}(4,i) = sigmaFLE(2,2);
        fleEstValues{j}(5,i) = sigmaFLE(2,3);
        fleEstValues{j}(6,i) = sigmaFLE(3,3);

        fleEigValues{j}(:,i) = sort(eig(fle));
        fleestRMS(i,j) = sqrt(trace(sigmaFLE));
        fleEstEigValues{j}(:,i) = sort(eig(sigmaFLE));
        frecovRMS(i,j) = sqrt(trace(frecov{j}));

        if(j == 1)
            %             % do the Kalman filter thing.
            %             %predict
            %             curestfrecov = prevfrecov;
            %             estP = prevP + Q*Q';
            %
            %             %update
            %             curfrecov = curestfrecov + K * (frevec - curestfrecov);
            %             curP = estP - K*estP;
            %             K = estP*inv(estP + R);
            %
            %             % set up for the next pass in the loop.
            %             prevfre = curfrecov;
            %             prevP = curP;
            %
            %             frenoisycov(i,:) = frevec;
            %             fresmoothcov(i,:) = curfrecov';
            %             fresmoothRMS(i) = sqrt(curfrecov(1) + curfrecov(4) + curfrecov(6));

            if((i > ((N/2)-5)) && (i < ((N/2)+5) ))
                frevec'
                [frecov{1}(1,1), frecov{1}(1,2), frecov{1}(1,3), frecov{1}(2,2), frecov{1}(2,3), frecov{1}(3,3)]
                %setxfrm.R
                %xfrm.rot
                %mrk
                %xfrmmrk
                %U
                %L
                V
            end
        end
    end
    fleRMS(i) = sqrt(trace(fle));
    meanFLE = meanFLE/nMrks;
    meanFLERMS(i) = sqrt(trace(meanFLE));
    meanFLEEigValues(:,i) = sort(eig(fle));

end

timeperframe = (cputime-starttime)/N;
fprintf( 'Time per Frame: %f ms or %f Hz\n', timeperframe*1000, 1/timeperframe);

%frenoisycov{1}(6000,:)

figure(1);
plot(1:N, fleRMS, 'k:');
hold on;
plot(1:N, frecovRMS(:,1), '-', 'Color', (1/3)*[1, 1, 1]);
plot(1:N, frecovRMS(:,2), ':', 'Color', (1/3)*[1, 1, 1]);
plot(1:N, frecovRMS(:,3), '-.', 'Color', (1/3)*[1, 1, 1]);
plot(1:N, frecovRMS(:,4), '--', 'Color', (1/3)*[1, 1, 1]);
hold off;
titlestring = sprintf('Case %d: Windowed FRE RMS (%d frames)', casestudy, winsize);
title(titlestring);
xlabel('Frame');
ylabel('RMS (mm)');
legend('True FLE RMS', 'Mrk #1', 'Mrk #2', 'Mrk #3', 'Mrk #4');

figure(2);
plot(1:N, fleRMS, 'k--', 'LineWidth',2);
hold on;
plot(1:N, fleestRMS(:,1), '-', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 2);
plot(1:N, fleestRMS(:,2), '-', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 1);
plot(1:N, fleestRMS(:,3), '-.', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 2);
plot(1:N, fleestRMS(:,4), '--', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 2);
plot(1:N, meanFLERMS, 'k', 'LineWidth', 2);
%plot(1:N, measfreRMS, '.b');
%plot(1:N, sqrt(2)*RRMS(:,1), 'c');
hold off;
titlestring1 = sprintf('Case %d: Estimated FLE RMS', casestudy);
titlestring = sprintf('(FRE computed with %d frames window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('RMS (mm)', 'fontsize', 18);
legend('True FLE RMS', 'Mrk #1', 'Mrk #2', 'Mrk #3', 'Mrk #4', ...
    'Mean of Markers 1-4', 'Location', 'SouthEast');
set(gca, 'fontsize', 14);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(3); % show the computation on Marker 1.
plot(1:N, fleRMS, 'k--', 'LineWidth', 2);
hold on;
plot(1:N, frecovRMS(:,1), '-', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 1);
plot(1:N, fresmoothRMS(:,1), '-.', 'Color', (1/3)*[1, 1, 1], 'LineWidth', 2);
plot(1:N, fleestRMS(:,1), 'k', 'LineWidth', 2);
%plot(1:N, meanFLERMS, 'c');
%plot(1:N, measfreRMS, '.b');
plot(1:N, sqrt(nMrks/(nMrks-2))*fresmoothRMS(:,1), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
titlestring1 = sprintf('Case %d: Windowed FRE and Estimated FLE for Marker #1', casestudy);
titlestring = sprintf('(FRE computed with %d frames window)', winsize);
title({titlestring1, titlestring}, 'fontsize', 18);
xlabel('Frame', 'fontsize', 18);
ylabel('RMS (mm)', 'fontsize', 18);
legend('True FLE RMS', 'Windowed FRE', 'Moving Average of FRE', 'Estimated FLE',...
    'Estimated FLE using Fitzpatrick 1998', 'Location', 'SouthEast');
set(gca, 'fontsize', 14);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(4); % plot the eigenvalues.
plot(1:N, fleEigValues{1}(1,:), 'r:');
hold on;
plot(1:N, fleEstEigValues{1}(1,:), 'r');
plot(1:N, fleEigValues{1}(2,:), 'g:');
plot(1:N, fleEstEigValues{1}(2,:), 'g');
plot(1:N, fleEigValues{1}(3,:), 'b:');
plot(1:N, fleEstEigValues{1}(3,:), 'b');
hold off;

figure(5);
plot(1:N, frecovRMS(:,1), 'b');
hold on;
plot(1:N, fresmoothRMS(:,1), 'g');
hold off;

figure(6);
for i = 1:6
    subplot(6,1,i);
    plot(1:N, frenoisycov{1}(:,i), 'b');
    hold on;
    plot(1:N, fresmoothcov{1}(:,i), 'g', 'LineWidth', 2);
    hold off;
end

figure(7);
plot(1:N, fleValues{1}(1,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(1,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.2]);
title('Estimated FLE Component: \sigma_{11}', 'fontsize', 20);
set(gca, 'fontsize', 18);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Sigma11-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Sigma11-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(8);
plot(1:N, fleValues{1}(2,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(2,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([-0.2 0.2]);
title('Estimated FLE Component: \sigma_{12}', 'fontsize', 20);
set(gca, 'fontsize', 18);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Sigma12-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Sigma12-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(9);
plot(1:N, fleValues{1}(3,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(3,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([-0.2 0.2]);
title('Estimated FLE Component: \sigma_{13}', 'fontsize', 20);
set(gca, 'fontsize', 18);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Sigma13-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Sigma13-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(10);
plot(1:N, fleValues{1}(4,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(4,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.2]);
title('Estimated FLE Component: \sigma_{22}', 'fontsize', 20);
set(gca, 'fontsize', 18);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Sigma22-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Sigma22-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(11);
plot(1:N, fleValues{1}(5,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(5,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([-0.2 0.2]);
title('Estimated FLE Component: \sigma_{23}', 'fontsize', 20);
set(gca, 'fontsize', 18);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Sigma23-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Sigma23-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);

figure(12);
plot(1:N, fleValues{1}(6,:), 'k--', 'LineWidth', 2);
hold on;
plot(1:N, fleEstValues{1}(6,:), '-', 'LineWidth', 2, 'Color', (2/3)*[1, 1, 1]);
hold off;
ylim([0 0.2]);
title('Estimated FLE Component: \sigma_{33}', 'fontsize', 20);
set(gca, 'fontsize', 18);
cd E:\docs\research\phd\publications\conference_papers\SPIE2009\fle-predict\images
filename = sprintf('Case%02d-FLE-estimated-Sigma33-Mrk01-%dframes.eps', casestudy, winsize);
print('-depsc2','-tiff','-r300', filename);
filename = sprintf('Case%02d-FLE-estimated-Sigma33-Mrk01-%dframes.jpg', casestudy, winsize);
print('-djpeg', filename);
