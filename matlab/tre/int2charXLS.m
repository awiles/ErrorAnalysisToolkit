function str = int2charXLS(n)
% returns the char representing column #n in excel-file

q = floor(n/26);
r = mod(n,26);

if r == 0
    r = 26;
    q = q - 1;
end

str = '';
if q ~= 0
    str = strcat(str,char(64+q));
end
if r ~= 0
    str = strcat(str,char(64+r));
end