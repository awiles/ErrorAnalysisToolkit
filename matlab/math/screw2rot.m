function R = screw2rot(screw)

kx = screw.axis(1);
ky = screw.axis(2);
kz = screw.axis(3);

cT = cos(screw.angle);
sT = sin(screw.angle);
vT = 1 - cT;

R = [ (kx*kx*vT + cT) (kx*ky*vT -kz*sT) (kx*kz*vT + ky*sT); ...
    (kx*ky *vT + kz*sT) (ky*ky*vT + cT) (ky*kz*vT -kx*sT); ...
    (kx*kz*vT - ky*sT) (ky*kz*vT + kx*sT) (kz*kz*vT + cT)];