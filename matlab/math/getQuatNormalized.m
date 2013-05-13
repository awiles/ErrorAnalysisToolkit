function q = getQuatNormalized(q0);
%#eml
qNorm = sqrt(sum(q0.^2));
if (q0(1) < 0)
    qNorm = -1*qNorm;
end

q = (1/qNorm) * q0;