function FLE = getQuadraticFLE(pos, rate, model)

nMrks = size(pos,1);
FLE = zeros(3,3,nMrks);

for i = 1:nMrks
    % check that the point is in the volume.
    switch(model)
        case 'cartesian'
            FLE(:,:,i) = diag([rate(1)*pos(i,1)^2 rate(2)*pos(i,2)^2 rate(3)*pos(i,3)^2]);
        case 'radial'
            r2 = pos(i,1)^2 + pos(i,2)^2 + pos(i,3)^2;
            FLE(:,:,i) = diag([rate(1)*r2 rate(2)*r2 rate(3)*r2]);
        otherwise
            error('Invalid model type given.');
    end
end