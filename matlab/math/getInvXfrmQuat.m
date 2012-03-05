 function xfrm = getInvXfrmQuat(xfrm0)

%invert the quaternion
xfrm.rot = getQuatNormalized(getInvQuat(xfrm0.rot));

%rotate the coordinate frame location and change direction.
xfrm.pos = (-1) * getRotPointQuat( xfrm.rot, xfrm0.pos);

