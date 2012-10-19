function [ color ] = getOrderedPlotColor( num )
%getOrderedPlotColor Return a unique color for each item.
%   For each number we return an appropriate color for the plots.
%   1 - 'Blue'
%   2 - 'Green'
%   3 - 'Red'
%   4 - 'Orange'
%   5 - 'Purple'
%   6 - 'Gray'
%   7 - 'DeepSkyBlue'
%   8 - 'Black'
%   9 - 'DeepPink'
%  10 - 'Gold'

switch(num)
    case 1
        color = rgb('Blue');
    case 2
        color = rgb('Green');
    case 3
        color = rgb('Red');
    case 4
        color = rgb('Orange');
    case 5
        color = rgb('Purple');
    case 6
        color = rgb('Gray');
    case 7
        color = rgb('DeepSkyBlue');
    case 8
        color = rgb('Black');
    case 9
        color = rgb('DeepPink');
    case 10
        color = rgb('Gold');
    otherwise
        warning('Number exceeded colors, returning blue');
        color = rgb('Blue');
end


end

