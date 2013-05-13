% script to transform random data...  eventually to be turned into a
% function.
% Given, mean and covariance matrix.

N = 1000;

mu = [0 0 0]';
Sigma = [0.1444         0   -0.1419;
        0    0.0911         0;
        -0.1419         0    0.1642];

[V, D] = eig(Sigma);

unitPoints = randn(3,N);
circPoints = getCirclePoints( 2, N );

%transform into new space.
randPoints = V * (D^.5) * unitPoints + repmat(mu, 1, N);

figure(1);
plotNDI(unitPoints', 'scatter');

subplot(2,2,1);
hold on;
plot( circPoints(1,:), circPoints(2,:), 'k');
axis equal;
hold off;

subplot(2,2,3);
hold on;
plot( circPoints(1,:), circPoints(2,:), 'k');
axis equal;
hold off;

subplot(2,2,1);
hold on;
plot( circPoints(1,:), circPoints(2,:), 'k');
axis equal;
hold off;


plotStochastic( randPoints, mu, Sigma );
return;
figure(2);
plotNDI(randPoints', 'scatter');

