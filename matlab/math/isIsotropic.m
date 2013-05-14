function [ bIso ] = isIsotropic( Sigma )
%isIsotropoic Determines if the matrix is isotropic FLE.
%   Sigma - FLE, TRE or other covariance matrix.
%   bIso    - boolean saying if it is isotropic.

eigv = eig(Sigma);

if( mean(eigv) == eigv(1) )
    bIso = 1;
else
    bIso = 0;
end

end

