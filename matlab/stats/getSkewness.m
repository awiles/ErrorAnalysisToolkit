function sk = getSkewness(data)

%skewness.
%m = getCentralMoments(data);
%s = m(3) / m(2)^(1.5);
[N, col] = size(data);
m = mean(data);
s = std(data);
m = m(ones(N,1),:);

m3= mean( (data - m).^3 );
sm2 = sqrt(mean( (data - m).^2 ));
sk = m3/(sm2^3);