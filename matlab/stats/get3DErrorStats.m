function errPos = get3DErrorStats(refPos, measPos, bFit)

if( size(refPos) ~= size(measPos) )
    errString = sprintf( 'The number of reference positions is different than the number of measured positions.\n'); 
    error(errString );
end

if( nargin > 2)
    if( bFit )
        [xfrm, refPos0, rmsError] = TransformRigidBody(refPos, measPos);
    end
else
    refPos0 = refPos;
end

errPos = measPos - refPos0;
errDist = (sum(errPos.^2,2)).^(.5);

%% mean errors

meanErr = mean([errPos errDist]);
stdErr  = std([errPos errDist]);
rmsErr  = sqrt(mean([errPos errDist].^2));
covErr = cov(errPos);

fprintf('\tRMS\t\t\tMean\t\tSt.Dev\t\n');
fprintf('x\t%f\t%f\t%f\n', rmsErr(1), meanErr(1), stdErr(1) );
fprintf('y\t%f\t%f\t%f\n', rmsErr(2), meanErr(2), stdErr(2) );
fprintf('z\t%f\t%f\t%f\n', rmsErr(3), meanErr(3), stdErr(3) );
fprintf('3D\t%f\t%f\t%f\n', rmsErr(4), meanErr(4), stdErr(4) );

xfrm