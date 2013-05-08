function [R, t, FRE] = point_register(X, Y, weight_sqr)

% This function performs point-based rigid body registration. It returns a
% rotation matrix, translation vector and the FRE. If provided with square
% of weights, it performs weighted version of registration and the FRE
% returned is the weighted FRE.

if nargin < 2
    error('At least two input arguments are required.');
end

[Ncoords Npoints] = size(X);
[Ncoords_Y Npoints_Y] = size(Y); 

if Ncoords ~= 3 | Ncoords_Y ~= 3
    error('Each argument must have exactly three rows.')
elseif (Ncoords ~= Ncoords_Y) | (Npoints ~= Npoints_Y)
    error('X and Y must have the same number of columns.');
elseif Npoints < 3
    error('X and Y must each have 3 or more columns.');
end

if nargin==2
    Xbar = mean(X,2);  % X centroid
    Ybar = mean(Y,2);  % Y centroid
    Xtilde = X-repmat(Xbar,1,Npoints); % X relative to centroid
    Ytilde = Y-repmat(Ybar,1,Npoints); % Y relative to centroid
    H = Xtilde*Ytilde';  % cross covariance matrix
    [U S V] = svd(H);    % U*S*V' = H
    R = V*diag([1, 1, det(V*U)])*U';
    t = Ybar - R*Xbar;
    if nargout == 3
        FREvect = R*X + repmat(t,1,Npoints) - Y;
        FRE = sqrt(mean(sum(FREvect.^2,1)));
    end
else
    weight_sqr = diag(diag(weight_sqr))'; % make sure it is a row vector
    w_norm = weight_sqr/sum(weight_sqr);
    Xbar = X*w_norm';
    Ybar = Y*w_norm';
    Xtilde = X-repmat(Xbar,1,Npoints); % X relative to centroid
    Ytilde = Y-repmat(Ybar,1,Npoints); % Y relative to centroid
    H = Xtilde*diag(weight_sqr)*Ytilde';
    [U S V] = svd(H);    % U*S*V' = H
    R = V*diag([1, 1, det(V*U)])*U';
    t = Ybar - R*Xbar;
    if nargout==3   % compute weighted FRE
        X1 = R*X + repmat(t,1,Npoints);
        FRE = sqrt(trace((X1 - Y)*diag(weight_sqr)*(X1-Y)'/Npoints));
    end
end

return