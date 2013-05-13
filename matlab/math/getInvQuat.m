function q = getInvQuat( q0 );

qNorm = sum(q0.^2);
q = (1/qNorm)*getConjQuat(q0);