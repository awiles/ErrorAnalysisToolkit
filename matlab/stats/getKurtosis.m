function k = getKurtosis(data)

%Kutosis excess.
%m = getCentralMoments(data);
%k = (m(4) / (m(2)^2)) - 3;
[N, col] = size(data);
m = mean(data);
s = std(data);
m = m(ones(N,1),:);

m4 = mean( (data - m).^4 );
m2 = mean( (data - m).^2 );
k = m4/(m2^2)-3;

