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
maxErr  = max([errPos errDist]);
covErr = cov(errPos);

fprintf('\tRMS\t\tMean\tSt.Dev\tMax\t\n');
fprintf('x\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n', rmsErr(1), meanErr(1), stdErr(1), maxErr(1) );
fprintf('y\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n', rmsErr(2), meanErr(2), stdErr(2), maxErr(2) );
fprintf('z\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n', rmsErr(3), meanErr(3), stdErr(3), maxErr(3) );
fprintf('3D\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n', rmsErr(4), meanErr(4), stdErr(4), maxErr(4) );

if( nargin > 2)
    if( bFit )
        fprintf('The data was registered for a best fit set of statistics.\n');
        fprintf('    Xfrm: %3.6f, %3.6f, %3.6f, %3.6f, %3.2f, %3.2f, %3.2f\n', ...
            xfrm.rot, xfrm.pos);
    end
end

plot3(refPos0(:,1), refPos0(:,2), refPos0(:,3), '.k');
hold on;
plot3(measPos(:,1), measPos(:,2), measPos(:,3), '.r');
hold off;