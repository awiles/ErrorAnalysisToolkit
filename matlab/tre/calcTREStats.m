function [RMS, SigmaCov, SigmaCovPA] = calcTREStats(FLE, mrkP)
%#eml
%**************************************************************************
% Determine Target Registration Error
%
%       Written by Andrew Wiles, March 18, 2005
%           - July 25, 2005, bug in distance from tip to principal axis fixed.
%           - May 23, 2007, changed to calcTRE and the original calcTRE was
%             renamed to be calcTRE.old.m
%           - March 18, 2008, added non-homogenous FLE models.
%
%**************************************************************************
%   The formula implemented here is based on the result provided by
%   Fitzpatrick et al., "Predicting Error in Rigid-Body Point-Base
%   Registration", IEEE Trans. on Medical Imaging, Vol. 17, No. 5, October
%   1998.  The actual formula is based on the special case provided by West
%   et al., "Designing Optically Tracked Instruments for Image-Guided
%   Surgery", IEEE Trans. on Medical Imaging, Vol. 23, No. 5, May 2004.
%
%   TRE is determined from the formula:
%           TRE^2 = FLE^2/N*(1 + 1/3 * (sum(k=1:3) d_k^2/f_k^2) )
%               where
%                   FLE is the fiducial localizer error
%                   N   is the number of markers
%                   d_k is the distance of the tip from the kth principal
%                       axis.*
%                   f_k is the rms distance of the markers from the kth
%                       principal axis.
%                       * note that principal axis is assumed to be the
%                       coordinate frame of the rigid body centered at the
%                       demeaned location of the markers.
%
%   INPUTS
%           FLE     the given fiducial localizer error
%           mrkP    the marker and tip positions, the last point listed
%                   will be assumed to be the tip location.
%**************************************************************************

%% check the FLE input and set up the error tensor.
if( ndims(FLE) == 2 )
    if( (max(size(FLE)) == 1) && (min(size(FLE)) == 1) )
        sigmaG = (FLE^2)/3 * eye(3);    %FLE converted to equal variances.
    elseif( (max(size(FLE)) == 3) && (min(size(FLE)) == 1) )
        sigmaG = diag(FLE);             %assuming the variances are input.
    elseif( (max(size(FLE)) == 3) && (min(size(FLE)) == 3) )
        sigmaG = FLE;                   %assuming full covariance matrix provided.
    else
        error('FLE matrix is incorrect size\n');
    end
    bHomogenous = 1;
    %fprintf('Homogenous FLE detected.\n'); %sigmaG
elseif( ndims(FLE) == 3 )
    sigmaG = FLE;
    bHomogenous = 0;
    %fprintf('Non-Homogenous FLE detected.\n'); %sigmaG
    %sigmaGSize = size(sigmaG)
    %mrkPSize = size(mrkP)
    if( size(sigmaG, 3) ~= (size(mrkP,1)-1) ) % subtract one because the target is on the bottom.
        error('Inhomogenous FLE matrix does not have the same number of covariance matrices as there are markers.\n');
    end
end

%% set up the fiducial markers and tip locations.
% extract the fiducials
x0   = mrkP(1:(end-1),:);
% extract the tip.
tip0 = mrkP(end,:);
% number of markers.
N = size(x0,1);

%demean the markers.
x = x0 - repmat(mean(x0),size(x0,1),1);
tip = tip0 - mean(x0);

%% transform the markers into principal axes.
[U L V] = svd(x,0);
Lambda = diag(L);
% get the lambda squared values
Lambda2 = (Lambda).^2;

x = x*V; %x = (V'*x')';
tip = tip*V; %tip = (V'*tip')';

% transform error tensor into principal axes.
if( bHomogenous )
    sigma = V' * sigmaG * V;
else
    sigma = zeros(size(sigmaG));
    for i = 1:N
        sigma(:,:,i) = V' * sigmaG(:,:,i) * V;
    end
end

%% TODO: insert code to print out covariance matrices.

%% compute the TRE.

if( bHomogenous )
    %calculate the RMS^2 using the three terms.
    term1 = 1/N*trace(sigma);

    term2 = 0;
    for i = 1:3
        for j = 1:3
            if( j ~= i)
                term2 = term2 + tip(i)^2*(Lambda2(i) * sigma(j,j) + Lambda2(j)*sigma(i,i))/...
                    (Lambda2(i) + Lambda2(j))^2;
            end
        end
    end

    term3 = 0;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                if( (j ~= i) && (k ~= i) && (k ~= j) )
                    term3 = term3 + (Lambda2(j)*sigma(i,k)*tip(i)*tip(k))...
                        /((Lambda2(i) + Lambda2(j))*(Lambda2(j) + Lambda2(k)));
                end
            end
        end
    end

    RMS2 = term1 + term2 + term3;

    RMS = sqrt(RMS2);

    % calculate the covariance matrix.
    SigmaCov = zeros(3);
    for i = 1:3
        for j = 1:3
            SigmaCov(i,j) = sigma(i,j)/N;
            rotTerm = 0;
            for k = 1:3
                for m = 1:3
                    if( (k ~= i) && (m ~= j) )
                        Z = Lambda2(k) * (k == m) * sigma(i,j)...
                            - Lambda2(k) * (k == j) * sigma(i,m)...
                            - Lambda2(i) * (i == m) * sigma(k,j)...
                            + Lambda2(i) * (i == j) * sigma(k,m);
                        rotTerm = rotTerm + (tip(k) * tip(m) * Z)...
                            / ((Lambda2(k) + Lambda2(i))*(Lambda2(m) + Lambda2(j)));
                    end
                end
            end
            SigmaCov(i,j) = SigmaCov(i,j) + rotTerm;
        end
    end

else %non-homogenous FLE.
    %[U0 L0 V0] = svd(x,0);
    U0 = U;
    L0 = L;
    V0 = eye(3);
    Lambda0 = diag(L0);
    % get the lambda squared values
    Lambda02 = (Lambda0).^2;

    %RMS first.
    term1 = 0;
    term2 = 0;
    term3 = 0;

    %term 1
    for a = 1:N
        for i = 1:3
            term1 = term1 + sigma(i,i,a);
        end
    end
    term1 = term1;

    %term 2
    for a = 1:N
        for i = 1:3
            for j = 1:3
                if ( i ~= j )
                    term2 = term2 +...
                        tip(j)*(Lambda0(j)*U0(a,j)*sigma(i,i,a) - Lambda0(i)*U0(a,i)*sigma(j,i,a))...
                        /(Lambda02(j) + Lambda02(i));
                end
            end
        end
    end
    term2 = term2;

    %term 3
    for a = 1:N
        for i = 1:3
            for j = 1:3
                for k = 1:3
                    if ( (i ~= j) && (i ~= k) )
                        term3 = term3 + tip(j)*tip(k)*(Lambda0(j)*Lambda0(k)*U0(a,j)*U0(a,k)*sigma(i,i,a) ...
                            - Lambda0(j)*Lambda0(i)*U0(a,j)*U0(a,i)*sigma(i,k,a) ...
                            - Lambda0(i)*Lambda0(k)*U0(a,i)*U0(a,k)*sigma(j,i,a) ...
                            + Lambda0(i)*Lambda0(i)*U0(a,i)*U0(a,i)*sigma(j,k,a)) ...
                            / ((Lambda02(j) + Lambda02(i))*(Lambda02(k) + Lambda02(i)));
                    end
                end
            end
        end
    end

    % add them up and take the sqrt.
    %     term1
    %     term2
    %     term3
    RMS2 = (1/N^2)*term1 + (2/N)*term2 + term3;
    RMS = sqrt(RMS2);

    % calculate the covariance matrix.
    sigma;
    SigmaCov = zeros(3);
    for i = 1:3
        for j = 1:3
            term1 = 0;
            term2a = 0;
            term2b = 0;
            term3 = 0;
            for a = 1:N
                term1 =  term1 + sigma(i,j,a);
                for k = 1:3
                    if( k ~= i)
                        term2a = term2a + (tip(k))*(Lambda0(k)*U0(a,k)*sigma(i,j,a) - Lambda0(i)*U0(a,i)*sigma(k,j,a))...
                            /(Lambda02(k) + Lambda02(i));
                        for m = 1:3
                            if( (m ~= j) )
                                term3 = term3 + tip(k)*tip(m)*(Lambda0(k)*Lambda0(m)*U0(a,k)*U0(a,m)*sigma(i,j,a) ...
                                    - Lambda0(k)*Lambda0(j)*U0(a,k)*U(a,j)*sigma(i,m,a)...
                                    - Lambda0(i)*Lambda0(m)*U0(a,i)*U(a,m)*sigma(k,j,a)...
                                    + Lambda0(i)*Lambda0(j)*U0(a,i)*U(a,j)*sigma(k,m,a))...
                                    / ((Lambda02(k) + Lambda02(i))*(Lambda02(m) + Lambda02(j)));
                            end
                        end
                    end
                    if( k ~= j)
                        term2b = term2b + tip(k)*(Lambda0(k)*U0(a,k)*sigma(j,i,a) - Lambda0(j)*U0(a,j)*sigma(k,i,a))...
                            /(Lambda02(k) + Lambda02(j));

                    end
                end
            end
            %term1
            %term2a
            %term2b
            %term3
            SigmaCov(i,j) = term1/(N^2) + (1/N)*(term2a + term2b) + term3;
        end
    end
end



%store the matrix in PA.
SigmaCovPA = SigmaCov;

%convert covariance matrix back into original space from principal axes.
SigmaCov = V * SigmaCov * V';
% rho = zeros(1,3);
% % check the covariance matrix for proper form.
% rho(1) = SigmaCov(1,2)/(sqrt(SigmaCov(1,1))*(sqrt(SigmaCov(2,2))));
% rho(2) = SigmaCov(2,3)/(sqrt(SigmaCov(2,2))*(sqrt(SigmaCov(3,3))));
% rho(3) = SigmaCov(1,3)/(sqrt(SigmaCov(1,1))*(sqrt(SigmaCov(3,3))));
% %if( sum(diag(SigmaCov) < 0) > 0 )
% if (min(diag(SigmaCov))<=0) 
%     warning('calcTRE::diagaonals on TRE covariance matrix are negative.');
%     fprintf('  SigmaCov   = |% 2.2f  % 2.2f  % 2.2f|\n', SigmaCov(1,:));
%     fprintf('               |% 2.2f  % 2.2f  % 2.2f|\n', SigmaCov(2,:));
%     fprintf('               |% 2.2f  % 2.2f  % 2.2f|\n\n', SigmaCov(3,:));
%     fprintf('  SigmaCovPA = |% 2.2f  % 2.2f  % 2.2f|\n', SigmaCovPA(1,:));
%     fprintf('               |% 2.2f  % 2.2f  % 2.2f|\n', SigmaCovPA(2,:));
%     fprintf('               |% 2.2f  % 2.2f  % 2.2f|\n\n', SigmaCovPA(3,:));
% end
% if (abs(max(rho)) > 1)
%     warning('calcTRE::correlations in TRE covariance matrix are not valid.');
%     fprintf('rho = %2.2f, %2.2f, %2.2f\n', rho);
% end
