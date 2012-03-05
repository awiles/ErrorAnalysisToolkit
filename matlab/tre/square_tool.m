%mrk = [0 25 0; 0 -25 0; -10 10 0; 10 -10 0];
mrk = [-75, 50, 0; 0, 50, 0; 75, 50, 0; -25, 0, 0;...
    25, 0, 0; -75, -50, 0; 0, -50, 0; 75, -50, 0 ];
tip = [0 0 0];
FLE = 0.15;
tip_dist = -250:10:250;

%********
% test 1
%********
for i=1:length(tip_dist)
    mrk0 = [mrk; tip + [0 tip_dist(i) 0] ];
    tre(i) = calcTRE_NDI(FLE, mrk0,0);
end

figure(2);
plot(tip_dist, tre);
title('TRE for tool tip along the y-axis');
xlabel('Tool tip location on y-axis (mm)');
ylabel('TRE (mm)');

print -dtiff TRE_y-axis
%********
% test 2
%********

[X, Y] = meshgrid(-250:10:250, -250:10:250);
for i = 1:size(X,1)
    for j = 1:size(X,2)
        mrk0 = [mrk; tip + [X(i,j), Y(i,j), 0]];
        tre_plane(i,j) = calcTRE_NDI(FLE, mrk0,0);
    end
end

figure(3);
surf(X,Y,tre_plane);

title('TRE for tool tip in the xy-plane');
xlabel('Tool tip location on x-axis (mm)');
ylabel('Tool tip location on y-axis (mm)');
zlabel('TRE (mm)');
print -dtiff TRE_plane

%********
% test 3
%********

angle = 0:1:359;
max_tip_dist = max(tip_dist);
for i=1:length(angle)
    mrk0 = [mrk; tip + (max_tip_dist*[0 cosd(angle(i)) sind(angle(i))]) ];
    tre_polar(i) = calcTRE_NDI(FLE, mrk0,0);
end

figure(4);
plot(angle, tre_polar);
title('TRE for tool tip rotated about the x-axis');
xlabel('Tool tip rotated about the x-axis (degrees)');
ylabel('TRE (mm)');