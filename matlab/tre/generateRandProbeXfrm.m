function probe = generateRandProbeXfrm(probe, maxAngle)

bValid = 0;

while(~bValid)
    % assume we get a valid orientation.
    bValid = 1;
    
    % define the transformation with a random orientation.
    xfrm.pos = [0 0 -1000];
    [xfrm.R, q] = getRandOrientation();
    
    
    % transform rigid body into test space.
    probe.Actual.mrk = (xfrm.R * probe.Rigid.mrk')' + repmat(xfrm.pos, 4, 1);
    probe.Actual.normals = (xfrm.R * probe.Rigid.normals')';
    probe.Actual.tip = (xfrm.R * probe.Rigid.tip')' + xfrm.pos;

    % check normals.
    for i = 1:size(probe.Actual.normals,1)
        ang = acosd(probe.Actual.normals(i,:)* [0 0 1]');
        if( (ang > maxAngle) || (ang < 0))
            %testResult = sprintf('Marker Normals Out of Range, theta = %3.2f',ang);
            bValid = 0;
            break;
        end
    end
end

probe.xfrm = xfrm;

