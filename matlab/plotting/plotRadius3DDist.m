function plotRadius3DDist(Sigma, t, lineprop)

step = t(2) - t(1);
t = t + (step/2)*ones(size(t));
g = radius3Ddist(t, Sigma);
plot(t, g, lineprop);