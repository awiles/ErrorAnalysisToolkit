function xfrm = readXfrmFile(filename,row)

inXfrms = csvread(filename);
nComponents = size(inXfrms,2);
if (nComponents ~= 7)
    error('File does not have 7 columns for the transormation data using quaternions.');
end

if(nargin < 2)
    row = 1;
end

xfrm.pos = inXfrms(row, 1:3);
xfrm.rot = inXfrms(row, 4:7);