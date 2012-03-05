function val = getPercentile(data, p)

if( (p<0) || (p>1)) 
    error('Invalid percentile value given.  Must be between 0 and 1.');
end
% assuming the data is a vector.
data = sort(data);
N = length(data);

kd = p*(N+1);
k = floor(kd);
d = kd - k;

if(k <=0) 
    val = data(1);
elseif (k >= N)
    val = data(N);
else
    val = data(k) + d*(data(k+1) - data(k));
end
