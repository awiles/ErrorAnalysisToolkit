function sigma = getFLEMatrix(parms, point)
%#eml
model = parms.modeltype;
switch model
    case 1 % radial model.
        vec = point - parms.origin;
        r = sqrt(vec * vec');
        rms = parms.rate * r + parms.baseline;
        rms2 = rms^2;
        sigma = weightMatrix(rms2, parms.weight);
    case 2 % random FLE.
        sigma = diag(rand(1,3) .* parms.range);
    case 3
        sigma = diag(rand(1,3) .* parms.range);
        R = getRandOrientation();
        sigma = R * sigma * R';
    otherwise
        error('Not a valid model given.');
end

