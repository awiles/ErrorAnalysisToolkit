function xfrm = getRigidXfrmSVD( refPos, measPos )
% computes the rigid body transform for the meas positions.
% based on the Schonemann's algorithm: Schonemann PH. "A generalized
% solution of the orthogonal Procrustes problem". Psychometrika, vol. 31,
% pp 1-10, 1966.  (Summarized in Fitzpatrick et al. TMI 17(5): 694-702, 1998).
%#eml
nRefMarkers = size(refPos, 1);
nMeasMarkers = size(measPos, 1);
if( nRefMarkers ~= nMeasMarkers)
    error('getRigidXfrm: Number of reference markers is not equal to the number of measured markers.\n');
end

%% find the mean positions for both the reference and the measured.
avgRefPos = mean(refPos);
avgMeasPos = mean(measPos);

%% "demean" the collection of points.
dmRefPos = refPos - repmat(avgRefPos, size(refPos,1), 1);
dmMeasPos = measPos - repmat(avgMeasPos, size(measPos,1), 1);

%% the SVD method.
% compute matrix YTX
YTX = dmMeasPos' * dmRefPos;
% get the SVD.
[A,D,B] = svd(YTX);
% compute rotation matrix.
R = B * diag([1,1,det(B*A)])* A';
% get the translation.
t = avgMeasPos - avgRefPos*R;
% fill structure and return.
xfrm.rot = R';
xfrm.pos = t;