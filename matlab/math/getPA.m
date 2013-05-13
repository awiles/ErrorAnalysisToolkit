function [V,D] = getPA(x0)

N = size(x0,1);
%%demean the markers.
x = x0 - repmat(mean(x0),N,1);

%get the moment of inertia tensor.
I0 = getmoi(x);             % see function below

% transform the markers into principal axes.
[V D] = eig(I0);
%x = (V'*x')';

function I0 = getmoi(x)

%find the inertia tensor.
I0 = zeros(3,3);

for i = 1:3
    for j=1:3
        if (i==j)
            for k=1:size(x,1)
                I0(i,j) = I0(i,j) + (sum(x(k,:).^2) - x(k,i)^2);
            end
        else
            for k=1:size(x,1)
                I0(i,j) = I0(i,j) - x(k,i)*x(k,j);
            end
        end
    end
end