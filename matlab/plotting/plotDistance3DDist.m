function plotDistance3DDist(S, t, lineformat, linecolor)

if(nargin < 3)
    lineformat = 'r-';
end

sxyz = sqrt(diag(S))
rho(1) = S(1,2)/sqrt(S(1,1)*S(2,2));
rho(2) = S(2,3)/sqrt(S(2,2)*S(3,3));
rho(3) = S(1,3)/sqrt(S(1,1)*S(3,3))

step = t(2) - t(1);
t = t + (step/2)*ones(size(t));
g = zeros(size(t));
for i=1:length(t)
    g(i) = edsvcd(t(i),sxyz, rho);
end
plot(t, g, lineformat,'Color', linecolor, 'LineWidth', 2);