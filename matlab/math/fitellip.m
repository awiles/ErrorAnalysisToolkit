% fitellip gives the 6 parameter vector of the algebraic circle fit
% to a(1)x^2 + a(2)xy + a(3)y^2 + a(4)x + a(5)y + a(6) = 0
% X & Y are lists of point coordinates and must be column vectors.
function a = fitellip(X,Y)

   % normalize data
   mx = mean(X);
   my = mean(Y);
   sx = (max(X)-min(X))/2;
   sy = (max(Y)-min(Y))/2;
   x = (X-mx)/sx;
   y = (Y-my)/sy;
   
   % Build design matrix
   D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

   % Build scatter matrix
   S = D'*D;

   % Build 6x6 constraint matrix
   C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

   % Solve eigensystem
   [gevec, geval] = eig(S,C);

   % Find the negative eigenvalue
   [NegR, NegC] = find(geval < 0 & ~isinf(geval));
   % sometimes the perfect data gives a zero eigenvalue.
   if(isempty(NegC))
       [NegR, NegC] = find(geval <= 0 & ~isinf(geval));
   end
   
   % Extract eigenvector corresponding to positive eigenvalue
   A = gevec(:,NegC);

   % unnormalize
   a = [
        A(1)*sy*sy,   ...
        A(2)*sx*sy,   ...
        A(3)*sx*sx,   ...
        -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
        -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
        A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
                - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
                + A(6)*sx*sx*sy*sy   ...
       ]';
