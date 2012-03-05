function p = getLinePlaneIntersection(p1, p2, p0, nrm)

% if( (size(p1) ~= [1,3]) || (size(p2) ~= [1,3]) || (size(p0) ~= [1,3]) || (size(nrm) ~= [1,3]) )
%     error('inputs are in the wrong dimension');
% end

u = (nrm * (p0 - p1)')/(nrm*(p2 - p1)');

if( (u < 0) || (u > 1) )
    % only return points between the segment ends.
    p = [nan nan nan];
else
    p = p1 + u * (p2 - p1);
end
