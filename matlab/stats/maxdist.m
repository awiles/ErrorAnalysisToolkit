function g = maxdist(t, rms)

g = zeros(size(t));

sigma = rms/sqrt(3);

for i = 1:length(t)
    g(i) = sqrt(2/pi) * (t(i)^2/sigma^3) * exp(-(t(i)^2)/(2*sigma^2));
end