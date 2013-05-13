function rot = getRotMatrixd(euler)

euler = (pi/180) * euler;

rot = getRotMatrix( euler );