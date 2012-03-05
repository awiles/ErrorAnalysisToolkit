function xfrm = getCombXfrm(xfrm1, xfrm2)

%% multiply the quaternions
xfrm.rot = getQuatMultiply( xfrm2.rot, xfrm1.rot);

%% rotate the origin of the one coordinate frame.
xfrm.pos = getRotPointQuat( xfrm2.rot, xfrm1.pos );
xfrm.pos = xfrm.pos + xfrm2.pos;
