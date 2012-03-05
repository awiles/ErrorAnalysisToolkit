function qq = getQuatMultiply( qL, qR )
% multiply two quaternions together (see Shoemake).
%% TODO:
% note NDI uses a different version which I think reduces the number of
% multiplcations (12 vs. 16).  
%#eml

w = 1; x = 2; y = 3; z = 4;

%% multiply the quats.
qq = zeros(1,4);
qq(w) = qL(w)*qR(w) - qL(x)*qR(x) - qL(y)*qR(y) - qL(z)*qR(z);
qq(x) = qL(w)*qR(x) + qL(x)*qR(w) + qL(y)*qR(z) - qL(z)*qR(y);
qq(y) = qL(w)*qR(y) + qL(y)*qR(w) + qL(z)*qR(x) - qL(x)*qR(z);
qq(z) = qL(w)*qR(z) + qL(z)*qR(w) + qL(x)*qR(y) - qL(y)*qR(x);

%% second method

fA = (qL(w) + qL(x)) * (qR(w) + qR(x));
fB = (qL(z) - qL(y)) * (qR(y) - qR(z));
fC = (qL(x) - qL(w)) * (qR(y) + qR(z));
fD = (qL(y) + qL(z)) * (qR(x) - qR(w));
fE = (qL(x) + qL(z)) * (qR(x) + qR(y));
fF = (qL(x) - qL(z)) * (qR(x) - qR(y));
fG = (qL(w) + qL(y)) * (qR(w) - qR(z));
fH = (qL(w) - qL(y)) * (qR(w) + qR(z));

qq(w) =  fB + (-fE - fF + fG + fH) / 2.0;
qq(x) =  fA - ( fE + fF + fG + fH) / 2.0;
qq(y) = -fC + ( fE - fF + fG - fH) / 2.0;
qq(z) = -fD + ( fE - fF - fG + fH) / 2.0;