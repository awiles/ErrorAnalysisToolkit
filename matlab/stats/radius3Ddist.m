function g = radius3Ddist(t, Sigma)
% t - distance locations.
% Sigma - covariance matrix.

g = zeros(size(t));
%diag(Sigma)
rho(1) = Sigma(1,2)/(sqrt(Sigma(1,1))*(sqrt(Sigma(2,2))));
rho(2) = Sigma(2,3)/(sqrt(Sigma(2,2))*(sqrt(Sigma(3,3))));
rho(3) = Sigma(1,3)/(sqrt(Sigma(1,1))*(sqrt(Sigma(3,3))));

for i = 1:length(t)
    g(i) = edsvcd(t(i), sqrt(diag(Sigma)), rho);
end