function p = chidist(x,k)

p = 0.5^(k/2)/gamma(k/2)*x.^(k/2 - 1) .* exp(-x/2);