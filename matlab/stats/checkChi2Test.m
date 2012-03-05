function result = checkChi2Test(chi2stat, dof)
% assume alpha = 0.05.
if( chi2stat < 0)
    error('Invalid chi2 statistic returned.');
end

switch(dof)
    case 6
        if( chi2stat < 12.59)
            result = 1;
        else
            result = 0;
        end
    case 9
        if( chi2stat < 16.92)
            result = 1;
        else
            result = 0;
        end
    otherwise
        error('Invalid dof returned.  Check the data.');
end