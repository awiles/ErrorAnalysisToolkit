function [rms, freVect] = getFRE(X,Y,xfrm)

if (size(X) ~= size(Y))
    error('X and Y must be the same size\n');
end

Xprime = zeros(size(X));

for i = 1:length(X)
    Xprime(i,:) = getXfrmPointQuat(xfrm, X(i,:));
end

freVect = Y - Xprime;
rms = sqrt(mean(sum(freVect.^2,2)));