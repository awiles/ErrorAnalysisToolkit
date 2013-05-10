function noise = getMeasNoise(fle, N, stream)
%#eml
noise = zeros(N,3);

if(ismatrix(fle))
    %[V,D] = eig(fle)
    [U,D,V] = svd(fle);
    noise = (randn(stream, [N, 3]))*(D.^0.5)*V';
    %fprintf('Identical FLE');
else
    for i = 1:N
        %[V,D] = eig(fle(:,:,i));
        [U,D,V] = svd(fle(:,:,i));
        noise(i,:) = (randn(stream, [1, 3]))*(D.^0.5)*V';
    end
end