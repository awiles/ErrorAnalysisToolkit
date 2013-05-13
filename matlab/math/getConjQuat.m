function q = getConjQuat(q0);
%#eml
q = zeros(1,4);
q(1) = q0(1);
q(2) = -q0(2);
q(3) = -q0(3);
q(4) = -q0(4);