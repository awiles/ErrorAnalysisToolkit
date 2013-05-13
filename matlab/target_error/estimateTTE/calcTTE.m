function [m,S,rms] = calcTTE(sigma,r,order)
x=1; y=2; z=3;

if(size(sigma,2) ~= size(sigma,1))
    error('calcTTE::The FLE matrix is not square.');
end

if(size(sigma,1) == 5)
    sensor = '5D';
    if( (r(x) ~= 0) || (r(y) ~=0))
        error('calcTTE::Tip locations for a 5D sensor must be along the z-axis.');
    end
elseif(size(sigma,1) == 6)
    sensor = '6D';
else
    error('calcTTE::The FLE matrix is the wrong size.  It needs to be 5x5 or 6x6.');
end

if(nargin < 3)
    order = 2;
end

S = zeros(3);
%% mean

switch(sensor)
    case '5D'
        switch(order)
            case 2 % all higher order terms are included.
                m =[0, 0, -2*(sigma(4,4) + sigma(5,5))*r(z)];
                %%cov_11
                S(1,1) = (-4 * sigma(4,4) * sigma(5,5) - 8 * sigma(4,5) ^ 2 ...
                    - 12 * sigma(5,5) ^ 2 + 4 * sigma(5,5)) * r(z) ^ 2 + 4 * r(z) * sigma(1,5) + sigma(1,1);
                %%cov_12
                S(1,2) = (12 * sigma(4,4) * sigma(4,5) + 12 * sigma(4,5) * sigma(5,5) ...
                    - 4 * sigma(4,5)) * r(z) ^ 2 + (-2 * sigma(1,4) + 2 * sigma(2,5)) * r(z) + sigma(1,2);
                S(2,1) = S(1,2);
                %%cov_13
                S(1,3) = 2 * r(z) * sigma(3,5) + sigma(1,3);
                S(3,1) = S(1,3);
                %%cov_22
                S(2,2) = (-12 * sigma(4,4) ^ 2 - 4 * sigma(4,4) * sigma(5,5) ...
                    - 8 * sigma(4,5) ^ 2 + 4 * sigma(4,4)) * r(z) ^ 2 - 4 * r(z) * sigma(2,4) + sigma(2,2);
                %%cov_23
                S(2,3) = -2 * sigma(3,4) * r(z) + sigma(2,3);
                S(3,2) = S(2,3);
                %%cov_33
                S(3,3) = (16 * sigma(4,5) ^ 2 + 8 * sigma(4,4) ^ 2 + 8 * sigma(5,5) ^ 2) * r(z) ^ 2 + sigma(3,3);
                %% RMS
                ms = (4 * sigma(4,4) + 4 * sigma(5,5)) * r(z) ^ 2 + (-4 * sigma(2,4) + 4 * sigma(1,5)) * r(z) + sigma(2,2) + sigma(1,1) + sigma(3,3);
                rms = sqrt(ms);
            case 1 % only the first order covariance terms are included.
                m =[0, 0, -2*(sigma(4,4) + sigma(5,5))*r(z)];
                %%cov_11
                S(1,1) = 4 * sigma(5,5) * r(z) ^ 2 + 4 * sigma(1,5) * r(z)  + sigma(1,1);
                %%cov_12
                S(1,2) = -4 * sigma(4,5) * r(z) ^ 2 + (-2 * sigma(1,4) + 2 * sigma(2,5)) * r(z) + sigma(1,2);
                S(2,1) = S(1,2);
                %%cov_13
                S(1,3) = 2 * sigma(3,5) * r(z)  + sigma(1,3);
                S(3,1) = S(1,3);
                %%cov_22
                S(2,2) = 4 * sigma(4,4) * r(z) ^ 2 - 4 * sigma(2,4) * r(z)  + sigma(2,2);
                %%cov_23
                S(2,3) = -2 * sigma(3,4) * r(z) + sigma(2,3);
                S(3,2) = S(2,3);
                %%cov_33
                S(3,3) = sigma(3,3);
                %% RMS
                ms = (4 * sigma(4,4) + 4 * sigma(5,5)) * r(z) ^ 2 + (-4 * sigma(2,4) + 4 * sigma(1,5)) * r(z) + sigma(2,2) + sigma(1,1) + sigma(3,3);
                rms = sqrt(ms);
            otherwise
                error('Invalid order given for the formulations.  Only order 1 (higher order terms ignored) and 2 (higher order terms included) are valid.');
        end
    case '6D'
        switch(order)
            case 2 % all higher order terms are included.
                m = [(-2 * sigma(5,5) - 2 * sigma(6,6)) * r(x) + 2 * sigma(4,5) * r(y) + 2 * sigma(4,6) * r(z),...
                    2 * sigma(4,5) * r(x) + (-2 * sigma(4,4) - 2 * sigma(6,6)) * r(y) + 2 * sigma(5,6) * r(z),...
                    2 * sigma(4,6) * r(x) + 2 * sigma(5,6) * r(y) + (-2 * sigma(4,4) - 2 * sigma(5,5)) * r(z)];
                %% cov_11.
                S(1,1) = (16 * sigma(5,6) ^ 2 + 8 * sigma(5,5) ^ 2 + 8 * sigma(6,6) ^ 2) * r(x) ^ 2 ...
                    + (r(y) * (-16 * sigma(5,5) * sigma(4,5) - 16 * sigma(4,6) * sigma(5,6))...
                    + r(z) * (-16 * sigma(6,6) * sigma(4,6) - 16 * sigma(4,5) * sigma(5,6))) * r(x)...
                    + r(y) ^ 2 * (4 * sigma(4,5) ^ 2 + 4 * sigma(4,4) * sigma(5,5) + 4 * sigma(6,6))...
                    + r(y) * (r(z) * (8 * sigma(4,5) * sigma(4,6) - 8 * sigma(5,6) + 8 * sigma(4,4) * sigma(5,6)) - 4 * sigma(1,6))...
                    + r(z) ^ 2 * (4 * sigma(4,4) * sigma(6,6) + 4 * sigma(4,6) ^ 2 + 4 * sigma(5,5)) ...
                    + 4 * r(z) * sigma(1,5) + sigma(1,1);
                %% cov_12
                S(1,2) = r(x) ^ 2 * (-8 * sigma(5,5) * sigma(4,5) - 8 * sigma(4,6) * sigma(5,6)) ...
                    + r(x) * (r(z) * (4 * sigma(4,4) * sigma(5,6) - 8 * sigma(6,6) * sigma(5,6)...
                    + 4 * sigma(4,5) * sigma(4,6) + 4 * sigma(5,6) - 8 * sigma(5,5) * sigma(5,6))...
                    + r(y) * (12 * sigma(4,5) ^ 2 + 4 * sigma(4,4) * sigma(5,5) + 8 * sigma(4,6) ^ 2 ...
                    - 4 * sigma(6,6) + 8 * sigma(6,6) ^ 2 + 8 * sigma(5,6) ^ 2) + 2 * sigma(1,6)) ...
                    + r(y) ^ 2 * (-8 * sigma(4,6) * sigma(5,6) - 8 * sigma(4,5) * sigma(4,4)) ...
                    + ((-8 * sigma(6,6) * sigma(4,6) - 8 * sigma(4,6) * sigma(4,4) + 4 * sigma(4,6) ...
                    + 4 * sigma(5,5) * sigma(4,6) + 4 * sigma(4,5) * sigma(5,6)) * r(z) - 2 * sigma(2,6)) * r(y) ...
                    + r(z) ^ 2 * (-4 * sigma(4,5) + 4 * sigma(4,6) * sigma(5,6) + 4 * sigma(6,6) * sigma(4,5)) ...
                    + (-2 * sigma(1,4) + 2 * sigma(2,5)) * r(z) + sigma(1,2);
                S(2,1) = S(1,2);
                %% cov_13
                S(1,3) = r(x) ^ 2 * (-8 * sigma(4,5) * sigma(5,6) - 8 * sigma(6,6) * sigma(4,6)) ...
                    + r(x) * (r(y) * (4 * sigma(4,4) * sigma(5,6) - 8 * sigma(6,6) * sigma(5,6) ...
                    + 4 * sigma(4,5) * sigma(4,6) + 4 * sigma(5,6) - 8 * sigma(5,5) * sigma(5,6)) ...
                    + r(z) * (4 * sigma(4,4) * sigma(6,6) + 8 * sigma(5,5) ^ 2 + 12 * sigma(4,6) ^ 2 ...
                    + 8 * sigma(5,6) ^ 2 + 8 * sigma(4,5) ^ 2 - 4 * sigma(5,5)) - 2 * sigma(1,5)) ...
                    + r(y) ^ 2 * (4 * sigma(4,5) * sigma(5,6) + 4 * sigma(5,5) * sigma(4,6) ...
                    - 4 * sigma(4,6)) + ((-8 * sigma(5,5) * sigma(4,5) - 8 * sigma(4,5) * sigma(4,4) ...
                    + 4 * sigma(4,6) * sigma(5,6) + 4 * sigma(4,5) + 4 * sigma(6,6) * sigma(4,5)) * r(z) ...
                    + 2 * sigma(1,4) - 2 * sigma(3,6)) * r(y) ...
                    + r(z) ^ 2 * (-8 * sigma(4,6) * sigma(4,4) - 8 * sigma(4,5) * sigma(5,6)) ...
                    + 2 * r(z) * sigma(3,5) + sigma(1,3);
                S(3,1) = S(1,3);
                %% cov_22
                S(2,2) = r(x) ^ 2 * (4 * sigma(4,5) ^ 2 + 4 * sigma(4,4) * sigma(5,5) + 4 * sigma(6,6)) ...
                    + (r(y) * (-16 * sigma(4,5) * sigma(4,4) - 16 * sigma(4,6) * sigma(5,6)) ...
                    + r(z) * (-8 * sigma(4,6) + 8 * sigma(5,5) * sigma(4,6) + 8 * sigma(4,5) * sigma(5,6)) + 4 * sigma(2,6)) * r(x) ...
                    + (16 * sigma(4,6) ^ 2 + 8 * sigma(4,4) ^ 2 + 8 * sigma(6,6) ^ 2) * r(y) ^ 2 ...
                    + r(y) * r(z) * (-16 * sigma(6,6) * sigma(5,6) - 16 * sigma(4,5) * sigma(4,6)) ...
                    + r(z) ^ 2 * (4 * sigma(4,4) + 4 * sigma(5,6) ^ 2 + 4 * sigma(5,5) * sigma(6,6)) ...
                    - 4 * sigma(2,4) * r(z) + sigma(2,2);
                %% cov_23
                S(2,3) = r(x) ^ 2 * (-4 * sigma(5,6) + 4 * sigma(4,5) * sigma(4,6) + 4 * sigma(4,4) * sigma(5,6)) ...
                    + ((-8 * sigma(6,6) * sigma(4,6) - 8 * sigma(4,6) * sigma(4,4) + 4 * sigma(4,6) ...
                    + 4 * sigma(5,5) * sigma(4,6) + 4 * sigma(4,5) * sigma(5,6)) * r(y) ...
                    + (-8 * sigma(5,5) * sigma(4,5) - 8 * sigma(4,5) * sigma(4,4) + 4 * sigma(4,6) * sigma(5,6) ...
                    + 4 * sigma(4,5) + 4 * sigma(6,6) * sigma(4,5)) * r(z) - 2 * sigma(2,5) + 2 * sigma(3,6)) * r(x) ...
                    + r(y) ^ 2 * (-8 * sigma(6,6) * sigma(5,6) - 8 * sigma(4,5) * sigma(4,6)) ...
                    + r(y) * (r(z) * (-4 * sigma(4,4) + 8 * sigma(4,6) ^ 2 + 4 * sigma(5,5) * sigma(6,6) ...
                    + 8 * sigma(4,5) ^ 2 + 8 * sigma(4,4) ^ 2 + 12 * sigma(5,6) ^ 2) + 2 * sigma(2,4)) ...
                    + r(z) ^ 2 * (-8 * sigma(5,5) * sigma(5,6) - 8 * sigma(4,5) * sigma(4,6)) ...
                    - 2 * r(z) * sigma(3,4) + sigma(2,3);
                S(3,2) = S(2,3);
                %% cov_33
                S(3,3) = r(x) ^ 2 * (4 * sigma(4,4) * sigma(6,6) + 4 * sigma(4,6) ^ 2 + 4 * sigma(5,5)) ...
                    + (r(y) * (8 * sigma(4,6) * sigma(5,6) - 8 * sigma(4,5) + 8 * sigma(6,6) * sigma(4,5)) ...
                    + r(z) * (-16 * sigma(4,5) * sigma(5,6) - 16 * sigma(4,6) * sigma(4,4)) - 4 * sigma(3,5)) * r(x) ...
                    + r(y) ^ 2 * (4 * sigma(4,4) + 4 * sigma(5,6) ^ 2 + 4 * sigma(5,5) * sigma(6,6)) ...
                    + r(y) * (r(z) * (-16 * sigma(5,5) * sigma(5,6) - 16 * sigma(4,5) * sigma(4,6)) + 4 * sigma(3,4)) ...
                    + (8 * sigma(5,5) ^ 2 + 8 * sigma(4,4) ^ 2 + 16 * sigma(4,5) ^ 2) * r(z) ^ 2 + sigma(3,3);
                %% rms
                ms = r(y) * (-4 * sigma(1,6) - 8 * sigma(5,6) * r(z) + 4 * sigma(3,4))...
                    + r(x) * (-4 * sigma(3,5) - 8 * sigma(4,6) * r(z) ...
                    + 4 * sigma(2,6) - 8 * sigma(4,5) * r(y))...
                    + (4 * sigma(5,5) + 4 * sigma(6,6)) * r(x) ^ 2 ...
                    + (4 * sigma(6,6) + 4 * sigma(4,4)) * r(y) ^ 2 ...
                    + (4 * sigma(5,5) + 4 * sigma(4,4)) * r(z) ^ 2 ...
                    + (-4 * sigma(2,4) + 4 * sigma(1,5)) * r(z) + sigma(1,1) + sigma(3,3) + sigma(2,2);
                rms = sqrt(ms);
            case 1 % only the first order covariance terms are included.
                m = [(-2 * sigma(5,5) - 2 * sigma(6,6)) * r(x) + 2 * sigma(4,5) * r(y) + 2 * sigma(4,6) * r(z),...
                    2 * sigma(4,5) * r(x) + (-2 * sigma(4,4) - 2 * sigma(6,6)) * r(y) + 2 * sigma(5,6) * r(z),...
                    2 * sigma(4,6) * r(x) + 2 * sigma(5,6) * r(y) + (-2 * sigma(4,4) - 2 * sigma(5,5)) * r(z)];
                %% cov_11.
                S(1,1) = 4 * sigma(6,6) * r(y) ^ 2 + r(y) * (-8 * sigma(5,6) * r(z) - 4 * sigma(1,6))...
                    +  4 * sigma(5,5)* r(z) ^ 2 + 4 * r(z) * sigma(1,5) + sigma(1,1);
                %% cov_12
                S(1,2) = r(x) * (4 * sigma(5,6) * r(z) - 4 * sigma(6,6) * r(y) + 2 * sigma(1,6)) ...
                    + (4 * sigma(4,6) * r(z) - 2 * sigma(2,6)) * r(y) - 4 * sigma(4,5) * r(z) ^ 2 ...
                    + (-2 * sigma(1,4) + 2 * sigma(2,5)) * r(z) + sigma(1,2);
                S(2,1) = S(1,2);
                %% cov_13
                S(1,3) = r(x) * (r(y) * 4 * sigma(5,6) - 4 * sigma(5,5)* r(z) - 2 * sigma(1,5)) ...
                    - 4 * sigma(4,6) * r(y) ^ 2 + (4 * sigma(4,5) * r(z) + 2 * sigma(1,4) - 2 * sigma(3,6)) * r(y) ...
                    + 2 * r(z) * sigma(3,5) + sigma(1,3);
                S(3,1) = S(1,3);
                %% cov_22
                S(2,2) = 4 * sigma(6,6) * r(x) ^ 2 + (-8 * sigma(4,6) * r(z) + 4 * sigma(2,6)) * r(x) ...
                    + 4 * sigma(4,4) * r(z) ^ 2 - 4 * sigma(2,4) * r(z) + sigma(2,2);
                %% cov_23
                S(2,3) = -4 * sigma(5,6) * r(x) ^ 2 + (4 * sigma(4,6) * r(y) + 4 * sigma(4,5) * r(z)...
                    - 2 * sigma(2,5) + 2 * sigma(3,6)) * r(x) + r(y) * (r(z) * -4 * sigma(4,4) + 2 * sigma(2,4)) ...
                    - 2 * r(z) * sigma(3,4) + sigma(2,3);
                S(3,2) = S(2,3);
                %% cov_33
                S(3,3) = 4 * sigma(5,5) * r(x) ^ 2 + (- 8 * sigma(4,5) * r(y) - 4 * sigma(3,5)) * r(x) ...
                    + 4 * sigma(4,4) * r(y) ^ 2 +  4 * sigma(3,4)* r(y) + sigma(3,3);
                %% rms
                ms = r(y) * (-4 * sigma(1,6) - 8 * sigma(5,6) * r(z) + 4 * sigma(3,4))...
                    + r(x) * (-4 * sigma(3,5) - 8 * sigma(4,6) * r(z) ...
                    + 4 * sigma(2,6) - 8 * sigma(4,5) * r(y))...
                    + (4 * sigma(5,5) + 4 * sigma(6,6)) * r(x) ^ 2 ...
                    + (4 * sigma(6,6) + 4 * sigma(4,4)) * r(y) ^ 2 ...
                    + (4 * sigma(5,5) + 4 * sigma(4,4)) * r(z) ^ 2 ...
                    + (-4 * sigma(2,4) + 4 * sigma(1,5)) * r(z) + sigma(1,1) + sigma(3,3) + sigma(2,2);
                rms = sqrt(ms);
            otherwise
                error('Invalid order given for the formulations.  Only order 1 (higher order terms ignored) and 2 (higher order terms included) are valid.');
        end
    otherwise
        error('calcTTE::Invalid sensor type.');
end

