FLE = 0.35;

x = [-175 0 90;
    0 18.61 -36.95;
    0 47.98 63.97;
    100 0 0];
pitch = 50;
R = [cos(pitch) 0 -sin(pitch);
    0 1 0;
    sin(pitch) 0 cos(pitch)];

xr = (R * x')'
return;

figure(1);
scatter3(x(:,1), x(:,2), x(:,3));
hold on;
scatter3(xr(:,1), xr(:,2), xr(:,3),'k');
hold off;
xlabel('x');
ylabel('y');
zlabel('z');

X = x - repmat(mean(x(1:end-1,:)),size(x,1),1)

[U,S,V] = svd(X(1:end-1,:),0)

Xp = [X] * V * inv(S)

figure(2)
scatter(x(1:3,3), x(1:3,1));
hold on;
scatter(xr(1:3,3), xr(1:3,1),'r');
scatter(x(4,3), x(4,1),'x');
scatter(xr(4,3), xr(4,1),'xr');
hold off;
xlabel('z');
ylabel('x');
axis([-10 110 -120 120]);
axis equal;
set(gca,'XDir','reverse'); 
set(gca,'YDir','reverse'); 

tip = [100, 0, 0];

TRE = calcTRE(FLE,[x;tip])



% [U0,S0,V0] = svd(X,0);
% 
% tail_1 = zeros(3); tip_1 = 10*eye(3);
% tail_2 = zeros(3); tip_2 = 10*V';
% 
% scatter3(X(:,1), X(:,2), X(:,3));
% hold on;
% quiver3(tail_1(:,1), tail_1(:,2), tail_1(:,3), tip_1(:,1), tip_1(:,2), tip_1(:,3), 'b');
% quiver3(tail_2(:,1), tail_2(:,2), tail_2(:,3), tip_2(:,1), tip_2(:,2), tip_2(:,3), 'g');
% hold off;
% 
% axis equal;

% x = [ 35.1791    -53.8979      9.0211;
%     -9.1391    -32.3543     -2.1145;
%     -13.8722     33.7664     -3.8523;
%     -65.6578     39.5779    -30.1377;
%     53.4901     12.9080     27.0834 ]
