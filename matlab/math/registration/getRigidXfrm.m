function xfrm = getRigidXfrm( refPos, measPos )
% computes the rigid body transform for the meas positions.
% based on the Horn algorithm: Horn, BKP. Closed-form solution of absolute
% orientation using unit quarternions. J. Opt. Soc. Am. A 4(4):629-642.
% April 1987.
%#eml
nRefMarkers = size(refPos, 1);
nMeasMarkers = size(measPos, 1);
if( nRefMarkers ~= nMeasMarkers)
    error('getRigidXfrm: Number of reference markers is not equal to the number of measured markers.\n');
end

% need to allocate for simulink.
xfrm.rot = [1 0 0 0];
xfrm.pos = [0 0 0];

%% find the mean positions for both the reference and the measured.
avgRefPos = mean(refPos);
avgMeasPos = mean(measPos);

%% "demean" the collection of points.
dmRefPos = refPos - repmat(avgRefPos, size(refPos,1), 1);
dmMeasPos = measPos - repmat(avgMeasPos, size(measPos,1), 1);

%% compute the matrix sums of products

% first the sums of products:
% is this just dmRefPos'*dmMeasPos to give the appropriate 3x3 matrix?
Sxx=0; Sxy=0; Syx=0; Sxz=0; Szx=0; Syy=0; Syz=0; Szy=0; Szz=0;
for i = 1:nRefMarkers
    Sxx = Sxx + dmRefPos(i,1)*dmMeasPos(i,1);
    Sxy = Sxy + dmRefPos(i,1)*dmMeasPos(i,2);
    Syx = Syx + dmRefPos(i,2)*dmMeasPos(i,1);
    Sxz = Sxz + dmRefPos(i,1)*dmMeasPos(i,3);
    Szx = Szx + dmRefPos(i,3)*dmMeasPos(i,1);
    Syy = Syy + dmRefPos(i,2)*dmMeasPos(i,2);
    Syz = Syz + dmRefPos(i,2)*dmMeasPos(i,3);
    Szy = Szy + dmRefPos(i,3)*dmMeasPos(i,2);
    Szz = Szz + dmRefPos(i,3)*dmMeasPos(i,3);
end

% test question above.
% M = dmRefPos' * dmMeasPos;
% Sxx = M(1,1); Sxy = M(1,2); Sxz = M(1,3);
% Syx = M(2,1); Syy = M(2,2); Syz = M(2,3);
% Szx = M(3,1); Szy = M(3,2); Szz = M(3,3);

% second build the N matrix.
N = zeros(4,4);
N(1,1) = Sxx + Syy + Szz;
N(1,2) = Syz - Szy; N(2,1) = N(1,2);
N(1,3) = Szx - Sxz; N(3,1) = N(1,3);
N(1,4) = Sxy - Syx; N(4,1) = N(1,4);
N(2,2) = Sxx - Syy - Szz;
N(2,3) = Sxy + Syx; N(3,2) = N(2,3);
N(2,4) = Szx + Sxz; N(4,2) = N(2,4);
N(3,3) = -Sxx + Syy - Szz;
N(3,4) = Syz + Szy; N(4,3) = N(3,4);
N(4,4) = -Sxx - Syy + Szz;

%% compute the eigenvalues and eigenvectors of the N matrix.
[V,D] = eig(N);

if( ~any(any(imag(V))) && ~any(any(imag(D))))
    [maxEigVal, maxEigId] = max(diag(real(D)));
    maxEigVec = real(V(:,maxEigId)');
    xfrm.rot = getQuatNormalized( maxEigVec );
else
    error('Complex eigenvalue decomposition.');
end

%% compute the best fit transformation.
rotAvgRefPos = getRotPointQuat( xfrm.rot, avgRefPos );
xfrm.pos = avgMeasPos - rotAvgRefPos;
