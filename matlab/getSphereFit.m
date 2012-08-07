function [p,r] = getSphereFit(X)
% this will fit a sphere to a point cloud.

[p1, r1] = getApproxSphereFit(X);

n = 0;
J = zeros(size(X,1), 4);
while(n < 100)
    for i=1:size(X,1)
        rmeas(i) = sqrt(sum((X(i,:) - p1).^2));
        d(i) = rmeas(i) - r1;
        %compute the Jacobian
        J(i,1) = -(X(i,1) - p1(1,1))/rmeas(i);
        J(i,2) = -(X(i,2) - p1(1,2))/rmeas(i);
        J(i,3) = -(X(i,3) - p1(1,3))/rmeas(i);
        J(i,4) = -1;
    end    
    
    %solve for the system paramters.
    u = inv(J'*J)*J' * -d';
    p1 = (p1' + u(1:3,1))';
    r1 = r1 + u(4,1);
    n=n+1;
end

p = p1;
r = r1;

function [p,r] = getApproxSphereFit(X)

A = zeros(size(X,1), 4);
b = zeros(size(X,1), 1);

for i = 1:size(X,1)
    A(i,1) = 2*X(i,1);
    A(i,2) = 2*X(i,2);
    A(i,3) = 2*X(i,3);
    A(i,4) = -1;
    b(i,1) = X(i,1)^2 + X(i,2)^2 + X(i,3)^2;
end

u = inv(A'*A)*A'*b;
p = u(1:3,1)';
r = sqrt( u(1,1)^2 + u(2,1)^2 + u(3,1)^2 -u(4,1));

