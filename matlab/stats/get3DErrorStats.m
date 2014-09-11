function errPos = get3DErrorStats(refPos, measPos, varargin)
%get3DErrorStats Compute the vector stats for a point cloud.
%   refPos  - N x 3 matrix of reference positions.
%   measPos - N x 3 matrix of the measured positions.

% default values for optional variables.
bFit = 0;
bPlot = 0;
verbose = 0;

if( nargin > 2 )
    nVarArgs = length(varargin);
    i = 1;
    while( i <= nVarArgs )
        if( strcmp(varargin{i}, 'BestFit'))
            bFit = 1;
        elseif( strcmp(varargin{i}, 'Plot'))
            bPlot = 1;
        elseif( strcmp(varargin{i}, 'Verbose'))
            verbose = 1;
        else
            warning('Unknown parameter: %s -- Ignoring it.', varargin{i})
        end
        i=i+1;
    end
end

if( size(refPos) ~= size(measPos) )
    errString = sprintf( 'The number of reference positions is different than the number of measured positions.\n');
    error(errString );
end

if( bFit )
    [xfrm, refPos0, rmsError] = TransformRigidBody(refPos, measPos);
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

fprintf('\tRMS\t\tMean\tSt.Dev\tMax\t\t2.5%%\t97.5%%\t95.0%%\n');
fprintf('x\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\tN/A\n', rmsErr(1), meanErr(1), stdErr(1), maxErr(1), getPercentile(errPos(:,1),0.025), getPercentile(errPos(:,1),0.975) );
fprintf('y\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\tN/A\n', rmsErr(2), meanErr(2), stdErr(2), maxErr(2), getPercentile(errPos(:,2),0.025), getPercentile(errPos(:,2),0.975) );
fprintf('z\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\tN/A\n', rmsErr(3), meanErr(3), stdErr(3), maxErr(3), getPercentile(errPos(:,3),0.025), getPercentile(errPos(:,3),0.975) );
fprintf('3D\t%3.2f\t%3.2f\t%3.2f\t%3.2f\tN/A   \tN/A   \t%3.2f\n', rmsErr(4), meanErr(4), stdErr(4), maxErr(4), getPercentile(errDist, 0.95) );
fprintf('\nCovariance Matrix:\n');
fprintf('% 3.2f\t% 3.2f\t% 3.2f\n', covErr);

if( bFit )
    fprintf('The data was registered for a best fit set of statistics.\n');
    fprintf('    Xfrm: %3.6f, %3.6f, %3.6f, %3.6f, %3.2f, %3.2f, %3.2f\n', ...
        xfrm.rot, xfrm.pos);
end

if( bPlot )
    plot3(refPos0(:,1), refPos0(:,2), refPos0(:,3), '.k');
    hold on;
    plot3(measPos(:,1), measPos(:,2), measPos(:,3), '.r');
    hold off;
    legend({'Reference','Measured'});
end
end