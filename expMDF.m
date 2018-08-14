% Write data to MDA dat file
% inputs: 
% 1- Data in table format 
% 2 - filename to write to *.dat
%
% usage: expMDF(dataTable, fileoutname)

% Author : Rahul Rajampeta
% Date: April 15, 2016

function expMDF(dataTable, fileoutname)
%%  create the fid to write
% fileoutname = 'C:\2014\MATLAB\Matlab_scripts\Mat_to_Dat\test.dat';
fid = fopen(fileoutname,'W');

% NOTE: dump empty characters 
% otherwise fseek wont move to desired position to write to
fwrite(fid, repmat(char(0),1,200000), 'char');

%%  pass the table data
% load('table_data.mat');
varNames = dataTable.Properties.VariableNames;
varUnits = getTableUnits(dataTable);
numChannels = length(varNames);
numSamples = height(dataTable);
rateSample = dataTable.time(2) - dataTable.time(1);

%% Define locations for each block
linkID = 0;                 % DO NOT CHANGE THIS
linkHD = 64;              % DO NOT CHANGE THIS
linkCC = 300; 
blocksize_cc = 63;      % for linear format of data 63 bytes
padd_cc = 2;
linkCN = linkCC + numChannels*(blocksize_cc+padd_cc);
blocksize_cn = 228;     % constant
padd_cn = 2;
linkCG = linkCN + numChannels*(blocksize_cn+padd_cn);
linkDG = linkCG + 50;   % 26+padding
linkDT = linkDG + 50;   %28 + padding

%% write ID block
writeIDBlock(fid, linkID);

%% write HD block
writeHDBlock(fid, linkHD, linkDG);

%% write DG block - assume only one DG
writeDGBlock(fid, linkDG, linkCG, linkDT);

%% write CG block - assume only one CG
writeCGBlock(fid, linkCG, linkCN, numChannels, numSamples);

%% write CC/CN blocks
for vars = 1 : numChannels
   writeCCBlock(fid, linkCC, varUnits{vars});
   nextCN = linkCN+blocksize_cn+padd_cn;
   writeCNBlock(fid, linkCN, linkCC, nextCN, varNames{vars}, vars, rateSample, vars == numChannels);
   
   linkCC = linkCC+blocksize_cc+padd_cc;
   linkCN = linkCN+blocksize_cn+padd_cn;
end

%%  and write DTblock
writeDTBlock(fid, linkDT, dataTable);

fclose(fid);
end

%% ID block def
% dont change this
function writeIDBlock(fid, offset)
fseek(fid,offset,'bof');
fwrite(fid,'MDF     3.00    TGT 15.0', 'char');
fwrite(fid,0,'uint16');
fwrite(fid,0,'uint16');
fwrite(fid,300,'uint16');
fwrite(fid,0,'uint16');
fwrite(fid,repmat(char(0),1,28),'char');
fwrite(fid,0,'uint16');
fwrite(fid,0,'uint16');
% ftell(fid)
end

%% HD block def
function writeHDBlock(fid, offset, linkDG)
blocksize = 164;

fseek(fid,offset,'bof');        % file pointer always to 64 for HD
fwrite(fid,char('HD'),'char');
fwrite(fid,blocksize,'uint16');
fwrite(fid,linkDG,'uint32');        % pointer to DG block
fwrite(fid,0,'uint32');                 % pointer to TX
fwrite(fid,0,'uint32');                 %pointer to PR
fwrite(fid,1,'uint16');                  % number of DG BLOCKS (assume only one to write to)
fwrite(fid,datestr(now,'dd:mm:yyyy'),'char');  % datestring
fwrite(fid,datestr(now,'HH:MM:SS'),'char');  % timestring
username = getenv('USERNAME');
fwrite(fid,username(1:min(length(username),31)), 'char');         % author
fseek(fid,32-length(username),0);     % fill the rest 32bytes with null 
fwrite(fid,repmat(char(' '),1,32),'char');          %organisation name
fwrite(fid,repmat(char(' '),1,32),'char');          % project name
fwrite(fid,repmat(char(' '),1,32),'char');          % vehicle/dyno name
end

%% DG block def
function writeDGBlock(fid, offset, linkCG, linkDT)
%assume only one DG block to write
blocksize = 28;

fseek(fid,offset,'bof');
fwrite(fid,char('DG'),'2*char');
fwrite(fid,blocksize,'uint16');        %block size
fwrite(fid,0,'uint32');     %skip next DG
fwrite(fid,linkCG,'uint32');    % CG link
fwrite(fid,0,'uint32');     % skip trigger block TR
fwrite(fid,linkDT,'uint32');        % Data DR link
fwrite(fid,1,'uint16');         % assume only one CG block
fwrite(fid,0,'uint16');         % assume zero recordIDs
fwrite(fid,0,'uint32');             %reserved
if blocksize ~= (ftell(fid)-offset)
    disp(['blocksize error DG']);
end
end

%% CG block def
function writeCGBlock(fid, offset, linkCN, numChannels, numSamples)
%assume only one CG block to write
blocksize = 26;

fseek(fid,offset,'bof');
fwrite(fid,char('CG'),'2*char');
fwrite(fid,blocksize,'uint16');    %block size
fwrite(fid,0,'uint32');     %skip next CG
fwrite(fid,linkCN,'uint32');        %link to CN block
fwrite(fid,0,'uint32');     %skip TX block
fwrite(fid,2,'uint16');     %recordID
fwrite(fid,uint16(numChannels),'uint16');     % number of channels
fwrite(fid,uint16(numChannels*4),'uint16');     % record size (number of channels*(single)=32bits=4bytes)
fwrite(fid,uint32(numSamples),'uint32');       %number of samples
if blocksize ~= (ftell(fid)-offset)
    disp(['blocksize error CG']);
end

end

%% CN block def
function writeCNBlock(fid, offset, linkCC, nextLinkCN, signalName, numVar, rateSample, lastVar)
if lastVar
    nextLinkCN = 0;
end
blocksize = 228;

fseek(fid,offset,'bof');
fwrite(fid,'CN','char');
fwrite(fid,blocksize,'uint16');    %block size
fwrite(fid,nextLinkCN,'uint32');        % link to nextCN
fwrite(fid,linkCC,'uint32');        % link to  CC block
fwrite(fid,0,'uint32');        % link to  CE block
fwrite(fid,0,'uint32');        % link to  CD block
fwrite(fid,0,'uint32');        % link to  TX block
if strcmp(signalName,'time')
    fwrite(fid,1,'uint16');     %channel type
else
    fwrite(fid,0,'uint16');
end
fwrite(fid,signalName(1:min(length(signalName),31)), 'char');
fseek(fid,32-length(signalName),0);     % fill the rest 32bytes with null 
fwrite(fid,'simulink output data','char');
fseek(fid,128-length('simulink output data'),0);     % fill the rest 128bytes with null 
fwrite(fid,32*(numVar-1),'uint16');    %start offset
fwrite(fid,32,'uint16');    %num of bits
fwrite(fid,2,'uint16');    % single format (4bytes)
fwrite(fid,0,'uint16');    % value range valid
fwrite(fid,0,'double');    % min signal value
fwrite(fid,0,'double');    % max signal value
fwrite(fid,rateSample,'double');    % sampling rate
fwrite(fid,0,'uint32');     %skip TX block
fwrite(fid,0,'uint32');     %skip TX block
fwrite(fid,0,'uint16');     % additional byte offset    
end

%% CC block def
function writeCCBlock(fid, offset, units)
% units modification - limit to 20 with trailing zeros
units = units(1:min(19,length(units)));
units = [units repmat(char(0),1,20-length(units))];

blocksize = 62;

fseek(fid,offset,'bof');
fwrite(fid,'CC','char');
fwrite(fid,blocksize,'uint16');    %block size
fwrite(fid,0,'uint16');    % value range valid
fwrite(fid,0,'double');    % min signal value
fwrite(fid,0,'double');    % max signal value
fwrite(fid,units,'char');    % units of the signal
fwrite(fid,0,'uint16');    % conversion type - linear
fwrite(fid,2,'uint16');    % number of parameters
fwrite(fid,0,'double');    % P1
fwrite(fid,1,'double');    % P2

end

%% Data block - write row-wise
function writeDTBlock(fid, offset, dataTable)
fseek(fid,offset,'bof');
% for row = 1:height(dataTable)
%     fwrite(fid,dataTable{row,:},'single');
% end
dataTableRow = reshape(dataTable{:,:}', 1, []);
fwrite(fid,dataTableRow, 'single');
end

%% DataTable units
function units = getTableUnits(dataTable)
units = dataTable.Properties.VariableUnits;

% if no units are set to either of the column
if isempty(units)
    % set empty units to first column, it forces empty units to all columns
    dataTable.Properties.VariableUnits{1} = '';
    units = dataTable.Properties.VariableUnits;
end
end
