% Aurora Comm
useTX = 0;
commPort = 'COM6';

% if we are using a virtual rom file, we need to open it.
romFileName = 'D:/tool-def-files/610029.rom';
[romFile, message] = fopen(romFileName);
if romFile == -1
    fprintf(message);
    error('File %s could not be opened.', romFileName);
end
% read in the file.
[romData,count] = fread(romFile,'uchar');
fclose(romFile);
if (count > 1024) %|| (count < 752)
    error('File %s could not be read properly.', romFileName);
end
% create the rom image.
romImg = zeros(2048,1);
for j=2:2:2048
    if j < 2*count
        bb = bitshift(romData(j/2),-4);
        romImg(j-1) = dec2hex(bb);
        cc = romData(j/2)-bitshift(bb,4);
        romImg(j) = dec2hex(cc);
    else
        romImg(j-1) = 48;
        romImg(j) = 48;
    end
end

% open up the serial port.
aurora = serial(commPort, 'BaudRate', 9600, 'terminator', 'CR', 'BytesAvailableFcnMode', 'terminator');
fopen(aurora);

%% send a serial break to ensure default set-up.
serialbreak(aurora)
[reply, count, msg] = fscanf(aurora);
if( count <= 0 )
    error('%s', msg)
elseif( ~strcmpi( reply(1:5), 'RESET') )
    fclose(aurora);
    error('Error: System did not receive the initial RESET');
end

%% reset the baud rate of the serial port.
% first change the aurora baud rate.
fprintf(aurora, 'COMM 60001');
reply = fscanf(aurora);
fprintf(['COMM 60001 Reply: ', reply]);
if( ~strcmpi(reply(1:4), 'OKAY') )
    fclose(aurora);
    error('Error: COMM 60001 failed');
end
% now change the system baud rate.
set(aurora, 'BaudRate', 921600, 'StopBits', 1, 'FlowControl', 'hardware');
pause(1);

%% init the system.
fprintf(aurora, 'INIT '); % note that the extra space is important here.
reply = fscanf(aurora);
fprintf(['INIT Reply: ', reply])
if( ~strcmpi(reply(1:4), 'OKAY') )
    fclose(aurora);
    error('Error: INIT command failed.');
end

%% set up the port handles.

% any port handles need freeing?
fprintf(aurora, 'PHSR 01');
phsr_reply = fscanf(aurora);
fprintf(['PHSR 01 Reply: ', phsr_reply]);
loc = 1;
nToFree = hex2dec(phsr_reply(1:2));
loc = loc + 2;
if(nToFree > 0)
    for i = 1:nToFree
        handle = phsr_reply(loc:loc+1);
        fprintf(aurora, 'PHF %s', handle);
        reply = fscanf(aurora);
        if( ~strcmpi(reply(1:4), 'OKAY') )
            fclose(aurora);
            error('Error: PHF command failed.');
        end
        loc = loc + 5;
    end
end

% write virtual rom file to port handle 01.
address = char('0000', '0040', '0080', '00C0', '0100', '0140', '0180', '01C0',...
    '0200', '0240', '0280', '02C0', '0300', '0340', '0380', '03C0');
ph = char('0A');
for j = 0:15
    out = num2str(romImg((j*128+1):((j+1)*128)),'%128c');
    send = sprintf('PVWR %2s%4s%s', ph, address(j+1,:), out);
    fprintf(aurora, send);
    fprintf(['Sending... ', send, '\n']);
    reply = fscanf(aurora);
    fprintf(['PVWR Reply: ', reply]);
    if( ~strcmpi(reply(1:4), 'OKAY') )
        fclose(aurora);
        error('Error: PVWR command failed.');
    end
end

% any port handles need to be initialized?
fprintf(aurora, 'PHSR 02');
phsr_reply = fscanf(aurora);
fprintf(['PHSR 02 Reply: ', phsr_reply]);
loc = 1;
nToInit = hex2dec(phsr_reply(1:2));
loc = loc + 2;
if(nToInit > 0)
    for i = 1:nToInit
        handle = phsr_reply(loc:loc+1);
        send = sprintf('PINIT %s', handle);
        fprintf(['Sending... ', send, '\n']);
        fprintf(aurora, send);
        reply = fscanf(aurora);
        if( ~strcmpi(reply(1:4), 'OKAY') )
            fclose(aurora);
            error('Error: PINIT command failed.');
        end
        loc = loc + 5;
    end
end

% any port handle need to be enabled?
fprintf(aurora, 'PHSR 03');
phsr_reply = fscanf(aurora);
fprintf(['PHSR 03 Reply: ', phsr_reply]);
loc = 1;
nToEnable = hex2dec(phsr_reply(1:2));
loc = loc + 2;
if(nToEnable > 0)
    for i = 1:nToEnable
        handle = phsr_reply(loc:loc+1);
        send = sprintf('PENA %sD', handle);
        fprintf(['Sending... ', send, '\n']);
        fprintf(aurora, send);
        reply = fscanf(aurora);
        if( ~strcmpi(reply(1:4), 'OKAY') )
            fclose(aurora);
            error('Error: PENA command failed.');
        end
        loc = loc + 5;
    end
end

% start tracking.
fprintf(aurora, 'TSTART ');
reply = fscanf(aurora);
fprintf(['TSTART Reply: ', reply]);
if( ~strcmpi(reply(1:4), 'OKAY') )
    fclose(aurora);
    error('Error: TSTART command failed.');
end


% get some data using TX..
for i = 1:50
    % tx
    fprintf(aurora, 'TX 0001');
    reply = fscanf(aurora);
    fprintf(['TX 0001 Reply: ', reply]);
end

% get some data using BX...
for i = 1:50
    % bx
    fprintf(aurora, 'BX 0001');
    fprintf('%s\n', aurora.TransferStatus);
    pause(1/40)
    if( aurora.BytesAvailable > 0 )
        reply = fread(aurora,aurora.BytesAvailable, 'char');
        hexreply = dec2hex(reply');
        bxreply = '';
        for j=1:length(hexreply)
            bxreply = [bxreply, hexreply(j,:)];
        end
        fprintf(['BX 0001 Reply: ', bxreply, '\n']);
    end
end



% stop tracking.
fprintf(aurora, 'TSTOP ');
reply = fscanf(aurora);
fprintf(['TSTOP Reply: ', reply]);
if( ~strcmpi(reply(1:4), 'OKAY') )
    fclose(aurora);
    error('Error: TSTOP command failed.');
end

fclose(aurora);
delete(aurora);
clear aurora;
