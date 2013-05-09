function A0 = xfrmToPAWrapper(x0)
%#eml

[U, L, V, x, A, Ainv] = xfrmToPA(x0);

% stack the matrices for now.
N = length(A);

A0 = zeros(N*6,6);
for i = 1:N
    A0(((i-1)*N+1):i*N,:) = A{i};
end