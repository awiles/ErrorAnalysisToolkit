function pos = getRotPointQuat( q, pos0);
%#eml

% %UcrossV = cross( q(2:4), pos0);
% UcrossV = [ q(3)*pos0(3) - q(4)*pos0(2),...
%             q(4)*pos0(1) - q(2)*pos0(3),...
%             q(2)*pos0(2) - q(3)*pos0(1) ];
% 
% pos = [ pos0(1) + 2.0 * (q(1)*UcrossV(1) + q(3)*UcrossV(3) - q(4)*UcrossV(2)),...
%         pos0(2) + 2.0 * (q(1)*UcrossV(2) + q(4)*UcrossV(1) - q(2)*UcrossV(3)),...
%         pos0(3) + 2.0 * (q(1)*UcrossV(3) + q(2)*UcrossV(2) - q(3)*UcrossV(1))];

r = [0 pos0];
rprime = getQuatMultiply(getQuatMultiply(q, r), getConjQuat(q));
pos = rprime(1,2:4);