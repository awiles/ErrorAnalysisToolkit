function probe = generateRandFiducialTargets(nMarkers, mrkRange, targetRange)

nValid = 0;
disThresh = 5;
angThresh = 2;
while(~nValid)
    probe.Rigid.mrk = (mrkRange/2)*(rand(nMarkers, 3) - 0.5);
    probe.Rigid.tip = (targetRange/2)*(rand(1, 3) - 0.5);

    % check for collinearity.
    mrkPairs = nchoosek(1:nMarkers,2);
    segPairs = nchoosek(1:size(mrkPairs,1),2);

    for i = 1:size(segPairs,1)
        vec1 = probe.Rigid.mrk(mrkPairs(segPairs(i,1),1),:)... 
            - probe.Rigid.mrk(mrkPairs(segPairs(i,1),2),:);
        d1 = sqrt(sum(vec1.^2,2));
        
        vec2 = probe.Rigid.mrk(mrkPairs(segPairs(i,2),1), :)... 
            - probe.Rigid.mrk(mrkPairs(segPairs(i,2),2), :);
        d2 = sqrt(sum(vec2.^2,2));
        theta = acosd( (vec1*vec2')/(sqrt(vec1*vec1')*sqrt(vec2*vec2')));
        if( d1 < disThresh )
            nValid = 0;
            warning('generateRandFiducialTargets:: distance threshold violated.');
            break;
        elseif( d2 < disThresh )
            nValid = 0;
            warning('generateRandFiducialTargets:: distance threshold violated.');
            break;
        elseif(theta < angThresh && theta > (180-angThresh)) 
            nValid = 0;
            warning('generateRandFiducialTargets:: angle threshold violated.');
            break;
        else
            nValid = 1;
        end
    end
end
