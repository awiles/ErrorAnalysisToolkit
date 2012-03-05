function euler = getEulerAngles(R)
% R is the rotation matrix.

Roll = atan2( R(2,1), R(1,1) );
CosRoll = cos( Roll );
SinRoll = sin( Roll );
euler.roll  = Roll;
euler.pitch = atan2( -R(3,1), CosRoll * R(1,1) + SinRoll * R(2,1) );
euler.yaw   = atan2( SinRoll * R(1,3)- CosRoll * R(2,3), ...
    -SinRoll * R(1,2) + CosRoll * R(2,2) );