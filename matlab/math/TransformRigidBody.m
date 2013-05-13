function [xfrm, xfrmPos, rmsError] = TransformRigidBody(refPos, measPos);

nRefMarkers = size(refPos, 1);
nMeasMarkers = size(measPos, 1);
if( nRefMarkers ~= nMeasMarkers)
    error('getRigidXfrm: Number of reference markers is not equal to the number of measured markers.\n');
end

%% compute the transform first.
xfrm = getRigidXfrm( refPos, measPos );

%% transform each of the reference positions into measured space and
%% compute the errors.
rmsError = 0;
nMarkers = 0;
for i = 1:size(refPos,1)
    xfrmPos(i,:) = getXfrmPointQuat(xfrm, refPos(i,:));
    rmsError = rmsError + sum((xfrmPos(i,:) - measPos(i,:)).^2);
    nMarkers = nMarkers + 1;
end

rmsError = sqrt(rmsError/nMarkers);