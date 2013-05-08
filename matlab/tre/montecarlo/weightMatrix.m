function sigma = weightMatrix(rms2, w)

rw2 = ((1/sum(w)) * w).^2;
sigma2 = (1/sum(rw2))*rms2;

sigma = diag(sigma2*rw2);