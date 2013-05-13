function rot = getRotMatrix(euler)

Rx = euler(1);
Ry = euler(2);
Rz = euler(3);

% from Craig 1989.
rotX = [ 1 0 0; 0 cos(Rx) -sin(Rx); 0 sin(Rx) cos(Rx)]; % yaw
rotY = [cos(Ry) 0 sin(Ry); 0 1 0; -sin(Ry) 0 cos(Ry)];  % pitch
rotZ = [cos(Rz) -sin(Rz) 0; sin(Rz) cos(Rz) 0; 0 0 1];  % roll

rot = rotZ*(rotY*rotX);