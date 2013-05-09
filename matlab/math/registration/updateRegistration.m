function [xfrm, fre] = updateRegistration(refmrk, measmrk, regmethod)
%% Written by Andrew D. Wiles, 2009-12-03.

nMrks = size(refmrk,1);

switch regmethod
    case 1
        % using SVD:
        xfrm = getRigidXfrmSVD( refmrk, measmrk );
        xfrmmrk = (xfrm.rot * refmrk')' + repmat(xfrm.pos, nMrks, 1);
    case 2
        % using Horn:
        %xfrm = getRigidXfrm( refmrk, measmrk );
        %xfrmmrk = (quat2rm(xfrm.rot) * refmrk')' + repmat(xfrm.pos, nMrks, 1);
    otherwise
        error('Invalid registration method specified.');
end
%compute the FRE.
fre = xfrmmrk - measmrk;