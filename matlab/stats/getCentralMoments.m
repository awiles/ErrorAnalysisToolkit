function m = getCentralMoments(data)

% for 5 central moments.
mu = mean(data);
N = length(data);
S = zeros(1,5);

for i = N
    d = (data(i) - mu);
    S(1) = S(1) + d;
    S(2) = S(2) + d^2;
    S(3) = S(3) + d^3;
    S(4) = S(4) + d^4;
    S(5) = S(5) + d^5;
end

%m(1) = 0;
%m(2) = (S(1)^2/N^2) + S(2)/N;
%m(3) = ((2*S(1)^3)/(N^3)) - (3*S(1)*S(2)/(N^2)) + S(3)/N;
%m(4) = -(3*(S(1)^4)/(N^4)) + ((6*(S(1)^2)*S(2))/(N^3)) - (4*S(1)*S(3))/(N^2) + (S(4))/(N);
m = (1/N) * S;