function g = normdist(t, mu, sigma)

g = zeros(size(t));

for i = 1:length(t)
    g(i) = (1/(sqrt(2*pi)*sigma)) * exp(-1*((t(i)- mu)^2)/(2*sigma));
end