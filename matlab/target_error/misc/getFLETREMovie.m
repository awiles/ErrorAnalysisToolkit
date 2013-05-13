%% create FLE/TRE Movie.
cd e:\temp
x = 2;
y = 3;
z = 1;


N = 25;
FLEMovie = moviein(N);
FLEFREMovie = moviein(N);
FLETREMovie = moviein(N);
TREMovie = moviein(N);

%% rigid body shape.
%test = 'neuro1';
test = 'West1';
switch(test)
    case 'neuro1'
        name = 'neuro1';
        FLE = eye(3);
        truept = [-95 -100 0; -120 0 0; -100 100 0; -70 105 0];
        target = [ -80 -10 0];
        axiswidth = 2.5;
        stack = 1;
    case 'West1'
        name = 'West1';
        FLE = 2*diag([0.11, 0.15, 0.33]);
        A = 71; B = 54; rho = 85;
        [truept, normals, target]=getWestToolDesign('e', A, B, rho);
        axiswidth = 2.5;
        stack = 0;
end

fig1axis = [min([truept(:,2);target(:,2)])-5 max([truept(:,2);target(:,2)])+5 ...
    min([truept(:,1);target(:,1)])-5 max([truept(:,1);target(:,1)])+5 ];

%% get the RMS.
[RMS, SigmaCov, SigmaCovPA] = calcTRE(FLE, [truept;target])

for i = 1:N    
    figure(1);
    set(gcf, 'position', [100 100 800 600]);
    if(stack)
        subplot(3,1,1);
    else
        subplot(1,3,1);
    end
    plot(truept(:,2), truept(:,1), 'og-', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    hold on;
    pt = truept + getMeasNoise(FLE, 4);
    plot(pt(:,2), pt(:,1), 'ob', 'MarkerSize',6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    hold off;
    titlestring = sprintf('Frame % 3d', i);
    title(titlestring);
    axis equal;
    axis(fig1axis);
    set(gca, 'YDir', 'reverse');
    
    %% plot the transformed points showing FRE.
    if(stack)
        subplot(3,1,2);
    else
        subplot(1,3,2);
    end
    
    xfrm = getRigidXfrm(truept,pt);
    xpt = zeros(size(pt));
    for j = 1:size(pt, 1)
        xpt(j,:) = getXfrmPointQuat(xfrm, truept(j,:));
    end
    xtarget = getXfrmPointQuat(xfrm, target);
    plot(truept(:,2), truept(:,1), 'og-', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    hold on;
    plot(target(:,2), target(:,1), 'oy', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    plot(xpt(:,2), xpt(:,1), 'om-', 'MarkerSize',6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    plot(xtarget(:,2), xtarget(:,1), 'oy', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    hold off;
    titlestring = sprintf('Frame % 3d', i);
    title(titlestring);
    axis equal;
    axis(fig1axis);
    set(gca, 'YDir', 'reverse');
    
    %% plot the transformed points showing FLE and FRE.
    if(stack)
        subplot(3,1,3);
    else
        subplot(1,3,3);
    end;
    
    plot(truept(:,2), truept(:,1), 'og-', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    hold on;
    plot(target(:,2), target(:,1), 'oy', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    plot(xpt(:,2), xpt(:,1), 'om-', 'MarkerSize',6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    plot(xtarget(:,2), xtarget(:,1), 'oy', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    plot(pt(:,2), pt(:,1), 'ob', 'MarkerSize',6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    hold off;
    titlestring = sprintf('Frame % 3d', i);
    title(titlestring);
    axis equal;
    axis(fig1axis);
    set(gca, 'YDir', 'reverse');
    
    FLETREMovie(i) = getframe(gcf);    
    
    figure(2); %FLE stuff.
    set(gcf, 'position', [100 100 800 800]);
    plotCovarianceRings(0.95, truept(1,:)', FLE, 'k', '-', 2, [x y z]);
    %plotCovarianceEllipse3RGB(0.95, 50, truept(1,:)', FLE, 2/3*ones(1,3), 1);
    hold on;
    plot3(truept(1,x), truept(1,y), truept(1,z), 'og-', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    plot3(pt(1,x), pt(1,y), pt(1,z), 'ob', 'MarkerSize',6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    hold off;
    axis([truept(1,x)-axiswidth truept(1,x)+axiswidth truept(1,y)-axiswidth truept(1,y)+axiswidth truept(1,z)-axiswidth truept(1,z)+axiswidth]);
    plotview = '3D';
    setNDIPlotViewProperties;
    FLEMovie(i) = getframe(gcf);
    
    figure(3); %FLE with FRE stuff.
    set(gcf, 'position', [100 100 800 800]);
    plotCovarianceRings(0.95, truept(1,:)', FLE, 'k', '-', 2, [x y z]);
    hold on;
    plot3(truept(1,x), truept(1,y), truept(1,z), 'og-', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    plot3(pt(1,x), pt(1,y), pt(1,z), 'ob', 'MarkerSize',6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    plot3(xpt(1,x), xpt(1,y), xpt(1,z), 'ob', 'MarkerSize',6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    quiver3(pt(1,x), pt(1,y), pt(1,z),xpt(1,x)-pt(1,x), xpt(1,y)-pt(1,y), xpt(1,z)-pt(1,z), 'k', 'LineWidth', 2, 'AutoScale', 'off');
    hold off;
    axis([truept(1,x)-axiswidth truept(1,x)+axiswidth truept(1,y)-axiswidth truept(1,y)+axiswidth truept(1,z)-axiswidth truept(1,z)+axiswidth]);
    plotview = '3D';
    setNDIPlotViewProperties;
    FLEFREMovie(i) = getframe(gcf);
    
    figure(4); %TRE stuff.
    set(gcf, 'position', [100 100 800 800]);
    plotCovarianceRings(0.95, target(1,:)', SigmaCov, 'k', '-', 2, [x y z]);
    hold on;
    plot3(target(1,x), target(1,y), target(1,z), 'oy', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    plot3(xtarget(1,x), xtarget(1,y), xtarget(1,z), 'oy', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    axis([target(1,x)-axiswidth target(1,x)+axiswidth target(1,y)-axiswidth target(1,y)+axiswidth target(1,z)-axiswidth target(1,z)+axiswidth]);
    hold off;
    plotview = '3D';
    setNDIPlotViewProperties;
    TREMovie(i) = getframe(gcf);
end

cd E:\docs\research\phd\experiments\FLEFRETRE-animations
mkdir(name);
cd(name);

movie2avi(FLETREMovie, 'FLETREMovie.avi', 'compression', 'none', 'quality', 100, 'fps', 5);
movie2avi(FLEMovie, 'FLEMovie.avi', 'compression', 'none', 'quality', 100, 'fps', 5);
movie2avi(FLEFREMovie, 'FLEFREMovie.avi', 'compression', 'none', 'quality', 100, 'fps', 5);
movie2avi(TREMovie, 'TREMovie.avi', 'compression', 'none', 'quality', 100, 'fps', 5);

figure(5);
plotCovarianceRings(0.95, [0 0 0]', FLE, 'b', '-', 2, [x y z]);
plotCovarianceEllipse3RGB(0.95, 50, [0 0 0]', FLE, 2/3*ones(1,3),1);
axis([-axiswidth axiswidth -axiswidth axiswidth -axiswidth axiswidth]);
print('-dpng', '-r300', 'FLE95.png');

figure(6);
plotCovarianceRings(0.95, [0 0 0]', SigmaCov, 'r', '-', 2, [x y z]);
plotCovarianceEllipse3RGB(0.95, 50, [0 0 0]', SigmaCov, 2/3*ones(1,3),1);
axis([-axiswidth axiswidth -axiswidth axiswidth -axiswidth axiswidth]);
print('-dpng', '-r300', 'TRE95.png');