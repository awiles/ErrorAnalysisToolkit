function plotSphere(center, radius, facecolor)
if(nargin<3)
    facecolor = 'r';
end
[XX, YY, ZZ] = ellipsoid(center(1),center(2),center(3),radius,radius,radius,50);
s = surface(XX,YY,ZZ);
set(s,'facecolor',facecolor,'edgecolor','none');
camlight left
lighting gouraud
alpha(0.5);
axis equal;