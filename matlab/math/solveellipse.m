%
% Given an ellipse in the form:
%       a(1)x^2 + a(2)xy + a(3)y^2 + a(4)x + a(5)y + a(6) = 0
% finds the standard form:
%       ((x-cx)/r1)^2 + ((y-cy)/r2)^2 = 1
% returning [r1 r2 cx cy theta]
%
function v = solveellipse(a)

        % get ellipse orientation
        theta = atan2(a(2),a(1)-a(3))/2;

        % get scaled major/minor axes
        ct = cos(theta);
        st = sin(theta);
        ap = a(1)*ct*ct + a(2)*ct*st + a(3)*st*st;
        cp = a(1)*st*st - a(2)*ct*st + a(3)*ct*ct;

        % get translations
        T = [[a(1) a(2)/2]' [a(2)/2 a(3)]'];
        t = -inv(2*T)*[a(4) a(5)]';
        cx = t(1);
        cy = t(2);

        % get scale factor
        val = t'*T*t;
        scale = 1 / (val- a(6));

        % get major/minor axis radii
        r1 = 1/sqrt(scale*ap);
        r2 = 1/sqrt(scale*cp);
        v = [r1 r2 cx cy theta]';