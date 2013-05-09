function fte = readFTECovariance(filename,row)

inCov = csvread(filename);
nComponents = size(inCov,2);
if (nComponents == 21)
    %6D
    dim = 6;
    fte = zeros(6);
elseif( size(inCov,2) == 15 )
    %5D
    dim = 5;
    fte = zeros(5);
else
    error('The width of the covariance file %s is not correct.', filename);
end

if(nargin < 2)
    row = 1;
end

nCnt = 1;
for i = 1:dim
    j = 1;
    while(j <= i)
        fte(i,j) = inCov(row,nCnt);
        fte(j,i) = fte(i,j);
        nCnt = nCnt+1;
        j= j+1;
    end
end
