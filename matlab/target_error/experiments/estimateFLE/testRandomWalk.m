close all;
clear all;

casestudy = 10;
body = 'west';
%body = 'west2';
%body = 'ta003-4';
cd E:\docs\research\phd\experiments\FLEPrediction\simulated-tool-paths

%% get the rigid body design
switch(body)
    case 'west'
        casename = 'WestTool';
        [refmrk, normals, tip] = getWestToolDesign('e', 71, 54, 85);
    case 'west2'
        casename = 'WestTool2';
        [refmrk, normals, tip] = getWestToolDesign('e', 71, 54, 85);
        refmrk = [refmrk; 0, 0, 5];
    case 'ta003-4'
        casename = 'TA003-4';
        refmrk = [0.0000, 0.0000, -0.0120;...
           -99.4066, -33.4878, 0.0417;...
           -153.6516, 0.2701, -0.0463;...
           -103.4595, 49.1463, 0.0287];
%         refmrk = [0.0000, 0.0000, 0.0;...
%             -99.4066, -33.4878, 0.0;...
%             -153.6516, 0.2701,  0.0;...
%             -103.4595, 49.1463, 0.0];
        tip = [200, 0, 0];
    otherwise
        error('Invalid body design.');
end
nMrks = size(refmrk,1);
%% pre-allocate memory and set up parameters.
N = 12000;
pos = zeros(N,3);
rot = zeros(N,3,3);
truemrks = cell(4);
%xfrm = cell(N);
for i = 1:nMrks
    truemrks{i} = zeros(N,3);
end
tooltip = zeros(N,3);
curpos = [0, 0, -1500];
curxfrm.pos = curpos;
curxfrm.R = eye(3);

step = [0, 0, 1];
stepsize = 0.1;
sigma = 0.01;
e_sigma = 0.1;

%% build a random walk path.
i = 1;
while i <= N
    % update the current pos.
    pos(i,:) = curpos;
    tooltip(i,:) = (curxfrm.R * tip')' + curxfrm.pos;
    mrk = (curxfrm.R * refmrk')' + repmat(curxfrm.pos, size(refmrk,1), 1);
    for j = 1:nMrks
        truemrks{j}(i,:) = mrk(j,:);
    end
    % save the current xfrm.
    xfrm{i} = curxfrm;

    % compute the next translation.
    lastpos = curpos;
    curpos = lastpos + stepsize*(step + sigma*randn(1,3));
    diff = curpos - lastpos;
    step = diff/sqrt(diff*diff');
    curxfrm.pos = curpos;
    % compute the next rotation.
    euler = e_sigma * randn(1,3);
    rot = getRotMatrixd(euler);
    curxfrm.R = rot*curxfrm.R;
    % increment the frame.
    i = i+1;
end

x = 2;
y = 3;
z = 1;
plot3(pos(:,x), pos(:,y), pos(:,z), 'b-');
hold on;
plot3(truemrks{1}(:,x), truemrks{1}(:,y), truemrks{1}(:,z), 'r-');
plot3(truemrks{2}(:,x), truemrks{2}(:,y), truemrks{2}(:,z), 'g-');
plot3(truemrks{3}(:,x), truemrks{3}(:,y), truemrks{3}(:,z), 'g-');
plot3(truemrks{4}(:,x), truemrks{4}(:,y), truemrks{4}(:,z), 'g-');
plot3(tooltip(:,x), tooltip(:,y), tooltip(:,z), 'y-');
for i = [1, 500:500:N]
    body = [truemrks{1}(i,:); truemrks{2}(i,:); truemrks{3}(i,:); truemrks{4}(i,:); truemrks{1}(i,:)];
    plot3(body(:,x), body(:,y), body(:,z), 'k-', 'linewidth', 2);
    shaft = [pos(i,:); tooltip(i,:)];
    plot3(shaft(:,x), shaft(:,y), shaft(:,z), 'k-', 'linewidth', 2);
    plot3(tooltip(i,x), tooltip(i,y), tooltip(i,z), 'k.', 'markersize', 12);
end
hold off;
axis equal;
xlabel('Y (mm)', 'fontsize', 18);
ylabel('Z (mm)', 'fontsize', 18);
zlabel('X (mm)', 'fontsize', 18);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
set(gca, 'fontsize', 18);

filename = sprintf('path-%02d', casestudy);
save(filename, 'xfrm');
filename = sprintf('%s-path-%02d', casename, casestudy);
print('-depsc2','-tiff','-r300', filename);
