function [mu, sigma, rms, t] = plotCumHist(data, t, lineprop)

N = length(data);
mu = mean(data);
sigma = std(data);
rms = sqrt(mean(data.^2));

if( isempty(t))
    minVal = min(data);
    maxVal = max(data);
    step = 0.01*(maxVal - minVal);
    t = minVal:step:maxVal;

else
    step = t(2) - t(1);    
end

d0 = histc(data, t);
d0cum = cumsum(d0);
t = t + (step/2)*ones(size(t));

plot(t, (1/N)*d0cum, lineprop);