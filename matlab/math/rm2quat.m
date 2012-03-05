function q = rm2quat(R)

w = 1; x = 2; y = 3; z = 4;

q(w) = 1 + R(1,1) + R(2,2) + R(3,3);
q(x) = 1 + R(1,1) - R(2,2) - R(3,3);
q(y) = 1 - R(1,1) + R(2,2) - R(3,3);
q(z) = 1 - R(1,1) - R(2,2) + R(3,3);

if( (q(w) > q(x)) && (q(w) > q(y)) && (q(w) > q(z)) )
    %q0 is the largest value.
    q(w) = sqrt(q(w)) / 2;
    if(q(w) == 0)
        return;
    end
    q(x) = (R(3,2) - R(2,3)) / (q(w) * 4);
    q(y) = (R(1,3) - R(3,1)) / (q(w) * 4);
    q(z) = (R(2,1) - R(1,2)) / (q(w) * 4);
elseif( (q(x) > q(w)) && (q(x) > q(y)) && (q(x) > q(z)) ) 
    %qx is the largest value.
    q(x) = sqrt(q(x)) / 2;
    if(q(x) == 0)
        return;
    end
    q(w) = (R(3,2) - R(2,3)) / (q(x) * 4);
    q(y) = (R(2,1) + R(1,2)) / (q(x) * 4);
    q(z) = (R(1,3) + R(3,1)) / (q(x) * 4);
elseif( (q(y) > q(w)) && (q(y) > q(x)) && (q(y) > q(z)) )
    %qy is the largest value.
    q(y) = sqrt(q(y)) / 2;
    if(q(y) == 0)
        return;
    end
    q(w) = (R(1,3) - R(3,1)) / (q(y) * 4);
    q(x) = (R(2,1) + R(1,2)) / (q(y) * 4);
    q(z) = (R(3,2) + R(2,3)) / (q(y) * 4);
else %( (q(z) > q(w)) && (q(z) > q(x)) && (q(z) > q(y)) )
    %qz is the largest value.
    q(z) = sqrt(q(z)) / 2;
    if(q(z) == 0)
        return;
    end
    q(w) = (R(2,1) - R(1,2)) / (q(z) * 4);
    q(x) = (R(1,3) + R(3,1)) / (q(z) * 4);
    q(y) = (R(3,2) + R(2,3)) / (q(z) * 4);
end

q = getQuatNormalized(q);