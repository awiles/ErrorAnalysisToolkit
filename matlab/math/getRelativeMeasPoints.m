function measPos = getRelativeMeasPoints(data, rotId, refId, toolId)
% returns the measured points in the frame of reference given by refId.

measPos = [];
for i = 1:size(data,1)
    j = rotId;
    for k = 1:size(data,3)
        pos = [];
        for m = 1:size(data{i,j,k,refId}.rawxfrm.pos,1)
            refXfrm.pos = data{i,j,k,refId}.rawxfrm.pos(m,:);
            refXfrm.rot = data{i,j,k,refId}.rawxfrm.quat(m,:);
            toolXfrm.pos = data{i,j,k,toolId}.rawxfrm.pos(m,:);
            toolXfrm.rot = data{i,j,k,toolId}.rawxfrm.quat(m,:);
            refXfrm = getInvXfrmQuat(refXfrm);
            xfrm = getCombXfrm(toolXfrm, refXfrm);
            pos = [pos; xfrm.pos];
        end
        measPos = [measPos; mean(pos)];
    end
end