function TRE = calcTREFitz(FLE, mrkP)

%**************************************************************************
% Determine Target Registration Error
% 
%       Written by Andrew Wiles, March 18, 2005
%           - July 25, 2005, bug in distance from tip to principal axis fixed.
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
%                   FLE is the fiducial localizer error FLE^2 = 3*stdev^2
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

% extract the fiducials
x0   = mrkP(1:(end-1),:);
% extract the tip.
tip0 = mrkP(end,:);
% number of markers.
N = size(x0,1);

%demean the markers.
x = x0 - repmat(mean(x0),size(x0,1),1);
tip = tip0 - mean(x0);

%get the moment of inertia tensor.
%I0 = getmoi(x);             % see function below

% transform the markers into principal axes.
%[V D] = eig(I0);
[U Lambda V] = svd(x);
x = x*V; %x = (V'*x')';
tip = tip*V; %tip = (V'*tip')';
[U0 Lambda0 V0] = svd(x);
%calculate rms distance of fiducials from each principal axis.
% bug fixed - July 25, 2005. The sum should not be here as the mean
% inherently adds to sum to the calculation.
% f(1) = sqrt(mean(sum(x(:,2).^2)+x(:,3).^2));
% f(2) = sqrt(mean(sum(x(:,1).^2)+x(:,3).^2));
% f(3) = sqrt(mean(sum(x(:,1).^2)+x(:,2).^2));
f(1) = sqrt(mean((x(:,2).^2)+x(:,3).^2));
f(2) = sqrt(mean((x(:,1).^2)+x(:,3).^2));
f(3) = sqrt(mean((x(:,1).^2)+x(:,2).^2));
% bug fixed - July 25, 2005.  We need the distance of the tip to
% principal axes not to the origin.
% d = sqrt(tip.^2); distance to local demeaned origin.
d(1) = sqrt(tip(:,2).^2 + tip(:,3).^2); %distance to first principal axis
d(2) = sqrt(tip(:,1).^2 + tip(:,3).^2); %distance to second principal axis
d(3) = sqrt(tip(:,1).^2 + tip(:,2).^2); %distance to third principal axis

% get ratios
% bug fixed - July 25, 2005.  Ratio needs to be squared.
% r = d./f;
r = (d./f).^2;

% use formula from Fitzpatrick.
TRE2 = (FLE^2/N)*(1+ 1/3*(sum(r)));

TRE = sqrt(TRE2);


function I0 = getmoi(x)

%find the inertia tensor.
I0 = zeros(3,3);

for i = 1:3
    for j=1:3
        if (i==j)
            for k=1:size(x,1)
                I0(i,j) = I0(i,j) + (sum(x(k,:).^2) - x(k,i)^2);
            end
        else
            for k=1:size(x,1)
                I0(i,j) = I0(i,j) - x(k,i)*x(k,j);
            end
        end
    end
end
