function xfrm = getNDIRelXfrm(refXfrm, toolXfrm)

%invRefXfrm = getInvXfrmQuat(refXfrm);
xfrm = getCombXfrm(toolXfrm, refXfrm);
%xfrm = getCombXfrm(invRefXfrm, toolXfrm);