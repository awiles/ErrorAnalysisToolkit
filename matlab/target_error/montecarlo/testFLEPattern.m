%% testFLEPattern.
clear all;
%close all;

[X Y Z] = meshgrid(-600:20:600, -600:20:600, -2000:20:0);

sigma_x = zeros(size(X));
sigma_y = zeros(size(X));
sigma_z = zeros(size(X));
sigma_rms = zeros(size(X));
%rate = [1.7778e-008, 1.7778e-008, 1.6000e-007];
rate = [5.7e-9, 5.7e-9, 5.1e-8];

for i = 1:size(X,1)
    for j = 1:size(X,2)
        for k = 1:size(X,3)
            pos = [X(i,j,k) Y(i,j,k) Z(i,j,k)];
            mat = getQuadraticFLE(pos, rate, 'radial');
            sigma_x(i,j,k) = sqrt(mat(1,1));
            sigma_y(i,j,k) = sqrt(mat(2,2));
            sigma_z(i,j,k) = sqrt(mat(3,3));
            sigma_rms(i,j,k) = sqrt(trace(mat));
        end
    end
end

%% plot the FLE RMS.
figure(1);
clf;
p1 = patch(isosurface(Y,Z,X,sigma_rms,0.1));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p1,'FaceColor', 0.75*[1,1,1], 'EdgeColor', 'none');
p2 = patch(isosurface(Y,Z,X,sigma_rms,0.2));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p2,'FaceColor', 0.5*[1,1,1], 'EdgeColor', 'none');
p3 = patch(isosurface(Y,Z,X,sigma_rms,0.3));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p3,'FaceColor', 0.25*[1,1,1], 'EdgeColor', 'none');
p4 = patch(isosurface(Y,Z,X,sigma_rms,0.4));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p4,'FaceColor', 0*[1,1,1], 'EdgeColor', 'none');
% daspect([1 1 1]);
% view(3); axis equal;
% camlight
% lighting gouraud
%camlight right
%lighting gouraud
alpha(0.9);
axis equal;
axis([-600, 600, -2000, 0, -600, 600]);
view([-37.5 30]);

%title('FLE RMS', 'fontsize', 20);
xlabel('y', 'fontsize', 18);
ylabel('z', 'fontsize', 18);
zlabel('x', 'fontsize', 18);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
legend('0.1 mm', '0.2 mm', '0.3 mm', '0.4 mm');
set(gca, 'fontsize', 18);
print('-depsc2','-tiff','-r300', 'FLEPatternRMS.eps');
print('-djpeg','-r300', 'FLEPatternRMS.jpg');

%% plot the FLE sigma_x.
figure(2);
clf;
p1 = patch(isosurface(Y,Z,X,sigma_x,0.02));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p1,'FaceColor', 0.75*[1,1,1], 'EdgeColor', 'none');
p2 = patch(isosurface(Y,Z,X,sigma_x,0.05));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p2,'FaceColor', 0.5*[1,1,1], 'EdgeColor', 'none');
p3 = patch(isosurface(Y,Z,X,sigma_x,0.1));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p3,'FaceColor', 0.25*[1,1,1], 'EdgeColor', 'none');
p4 = patch(isosurface(Y,Z,X,sigma_x,0.15));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p4,'FaceColor', 0*[1,1,1], 'EdgeColor', 'none');
% daspect([1 1 1]);
% view(3); axis equal;
% camlight
% lighting gouraud
%camlight right
%lighting gouraud
alpha(0.9);
axis equal;
axis([-600, 600, -2000, 0, -600, 600]);
view([-37.5 30]);

%title('\sigma_{x}', 'fontsize', 20);
xlabel('y', 'fontsize', 18);
ylabel('z', 'fontsize', 18);
zlabel('x', 'fontsize', 18);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
legend('0.02 mm', '0.05 mm', '0.10 mm', '0.15 mm');
set(gca, 'fontsize', 18);
print('-depsc2','-tiff','-r300', 'FLEPatternSigma11.eps');
print('-djpeg','-r300', 'FLEPatternSigma11.jpg');

%% plot the FLE sigma_y.
figure(3);
clf;
p1 = patch(isosurface(Y,Z,X,sigma_y,0.02));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p1,'FaceColor', 0.75*[1,1,1], 'EdgeColor', 'none');
p2 = patch(isosurface(Y,Z,X,sigma_y,0.05));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p2,'FaceColor', 0.5*[1,1,1], 'EdgeColor', 'none');
p3 = patch(isosurface(Y,Z,X,sigma_y,0.1));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p3,'FaceColor', 0.25*[1,1,1], 'EdgeColor', 'none');
p4 = patch(isosurface(Y,Z,X,sigma_y,0.15));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p4,'FaceColor', 0*[1,1,1], 'EdgeColor', 'none');
% daspect([1 1 1]);
% view(3); axis equal;
% camlight
% lighting gouraud
%camlight right
%lighting gouraud
alpha(0.9);
axis equal;
axis([-600, 600, -2000, 0, -600, 600]);
view([-37.5 30]);

%title('FLE \sigma_{y}', 'fontsize', 20);
xlabel('y', 'fontsize', 18);
ylabel('z', 'fontsize', 18);
zlabel('x', 'fontsize', 18);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
legend('0.02 mm', '0.05 mm', '0.10 mm', '0.15 mm');
set(gca, 'fontsize', 18);
print('-depsc2','-tiff','-r300', 'FLEPatternSigma22.eps');
print('-djpeg','-r300', 'FLEPatternSigma22.jpg');

%% plot the FLE sigma_z.
figure(4);
clf;
p1 = patch(isosurface(Y,Z,X,sigma_z,0.1));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p1,'FaceColor', 0.75*[1,1,1], 'EdgeColor', 'none');
p2 = patch(isosurface(Y,Z,X,sigma_z,0.2));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p2,'FaceColor', 0.5*[1,1,1], 'EdgeColor', 'none');
p3 = patch(isosurface(Y,Z,X,sigma_z,0.3));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p3,'FaceColor', 0.25*[1,1,1], 'EdgeColor', 'none');
p4 = patch(isosurface(Y,Z,X,sigma_z,0.4));
%isonormals(X,Y,Z,sigma_rms,p,'negate');
set(p4,'FaceColor', 0*[1,1,1], 'EdgeColor', 'none');
% daspect([1 1 1]);
% view(3); axis equal;
% camlight
% lighting gouraud
%camlight right
%lighting gouraud
alpha(0.9);
axis equal;
axis([-600, 600, -2000, 0, -600, 600]);
view([-37.5 30]);

%title('FLE \sigma_{z}', 'fontsize', 20);
xlabel('y', 'fontsize', 18);
ylabel('z', 'fontsize', 18);
zlabel('x', 'fontsize', 18);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
legend('0.1 mm', '0.2 mm', '0.3 mm', '0.4 mm');
set(gca, 'fontsize', 18);
print('-depsc2','-tiff','-r300', 'FLEPatternSigma33.eps');
print('-djpeg','-r300', 'FLEPatternSigma33.jpg');
