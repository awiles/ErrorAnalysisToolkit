function data = readFromXLS(filename)

% get number of lines in the file
fid = fopen(filename,'rt');
nLines = 0;
while (fgets(fid) ~= -1),
    nLines = nLines + 1;
end
fclose(fid);

% get number of tools
nTools = xlsread(filename,'A2:A2');

% initialization
data = cell(1,nTools); % data contains cells, each cell represents a tool
curPos = 1;

% parsing the data
for i = 1:nTools
    curPos = curPos + 5;
    range = strcat(int2charXLS(curPos),'2:',int2charXLS(curPos),'2');
    nMarkers = xlsread(filename,range);
    data{i} = cell(1,nMarkers);
    for j = 1:nMarkers
        curPos = curPos + 2;
        range = strcat(int2charXLS(curPos),'2:',int2charXLS(curPos+2),int2str(nLines));
        data{i}{j} = xlsread(filename,range);
        curPos = curPos + 2;
    end
end
