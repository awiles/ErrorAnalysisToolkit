function [chistat, chi2prob, dof] = WishartHypTest(mu0, Sigma0, data, test, alpha)

%% set the data up.
% get the dimensions of the data.
[N,p] = size(data);

% convert the data into the identity form.
% 1. demean the data.
datatemp = data - repmat(mu0, size(data,1), 1);
% 2. transform the data into the identity space.
SigmaInvSqrt = inv(sqrtm(Sigma0));
data1 = (SigmaInvSqrt*datatemp')';

% compute the observational statistics of the data.
mu1 = mean(data);
Sigma1 = cov(data);

V = data1'*data1;
%VoverN = V/N;

%% hypothesis test.

switch(test)
    case {'covariance'}
        % H0: Sigma = I
        % compute the lambda statistic (lambda^* in Srivastava)
        % lambda = (exp(1)/n)^(0.5*p*n)* det(V)^(0.5*n)*exp(-0.5*trace(V)); 
        % chistat = -2*log(lambda);
        n = N-1;
        chistat = -p*n + p*n*log(n) - n * log(det(V)) + trace(V);
        dof = 0.5*p*(p+1);
    case{'meanandcov'}
        % H0: mu = 0, Sigma = I
        % lambda = (exp(1)/N)^(0.5*p*N)* det(V)^(0.5*N)*exp(-0.5*trace(V + N*mu1'*mu1));
        % chistat = -2*log(lambda);
        chistat = -p*N + p*N*log(N) - N * log(det(V)) + trace(V + N * mu1'*mu1);
        dof = 0.5*p*(p+1) + p;
    otherwise
        error('Invalid Hypothesis Test Type Provided');
end

chi2prob = chidist(chistat,dof);

% if(nargin > 4)
%     switch(alpha)
%         case 0.10
%             
%         case 0.05
%         case 0.01
%         otherwise
%             error('Invalid alpha value provided.');
%     end
% end