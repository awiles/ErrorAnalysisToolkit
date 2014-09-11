function [rms, freVect] = getFRE(X,Y,xfrm, varargin)

if (size(X) ~= size(Y))
    error('X and Y must be the same size\n');
end

verbose = 0;

if( nargin > 3 )
    nVarArgs = length(varargin);
    i = 1;
    while( i <= nVarArgs )
        if( strcmp(varargin{i}, 'verbose') || strcmp(varargin{i}, 'Verbose'))
            verbose = 1;
        else
            warning('Unknown parameter: %s -- Ignoring it.', varargin{i})
        end
        i=i+1;
    end
end

Xprime = zeros(size(X));

for i = 1:length(X)
    Xprime(i,:) = getXfrmPointQuat(xfrm, X(i,:));
end

freVect = Y - Xprime;
rms = sqrt(mean(sum(freVect.^2,2)));

if(verbose)
    fprintf('Registration FRE: %f\n', rms);
end