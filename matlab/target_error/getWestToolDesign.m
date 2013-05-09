function [mrk, normals, tip]=getWestToolDesign(design, A, B, rho)
%#eml
switch(design)
    case {'d'}
        % West design in Fig. 2(d).
        mrk = [ 0 B/2 0;
            -A/2 0 0;
            0 -B/2 0;
            A/2 0 0 ];
        normals = repmat([ 0 0 1], 4, 1);
    case{'e'}
        % West design in Fig. 2(e).
        mrk = [ -A/2 B/2 0;
            -A/2 -B/2 0;
            A/2 -B/2 0;
            A/2 B/2 0 ];
        normals = repmat([ 0 0 1], 4, 1);
    otherwise
        %warning('Invalid tool design case: use 2(e).');
        % West design in Fig. 2(e).
        mrk = [ -A/2 B/2 0;
            -A/2 -B/2 0;
            A/2 -B/2 0;
            A/2 B/2 0 ];
        normals = repmat([ 0 0 1], 4, 1);
end

% tip is located at
tip = [ rho 0 0 ];