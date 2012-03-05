function pos = getXfrmPointQuat(xfrm, pos0)

pos = getRotPointQuat( xfrm.rot, pos0) + xfrm.pos;