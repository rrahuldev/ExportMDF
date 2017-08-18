function mat2dat(varargin)
% Uses variables stored in timeseries format and
% converts to -dat file to read it in MDA
% saves the file name as *_mda.dat in the same location 
disp('  ');
%% check matlab version : currently supports 8.2 and higher
if verLessThan('matlab', '8.2')
    disp(repmat(char('-'), 1, 90));
    disp(' Currently, this script supports only version 8.2 (R2013b) and higher ');
    disp(repmat(char('-'), 1, 90));
    disp('  ');
    return;
end

%% check invoke method of the script
switch nargin
    case 0
        % Select mat file with variables stored as timeseries
        % UI for file selection
        [FileName,PathName,~] = uigetfile(fullfile(pwd,'*.mat'), 'Select the MAT file(s). Use Ctrl key for selecting multiple files','MultiSelect','on');

        % check one or multiple files selected
        if isa(FileName, 'char')
            % only one file selected
            % convert to cell format
            FileName = {FileName};
        elseif ~isa(FileName, 'cell')
            % no file selected - end the function call 
            return
        end
        
        % load each file and create dat file for each
        for f = 1: length(FileName)
            disp('');disp(['Loading file -- "', FileName{f},'"']);
            raw = load(strcat(PathName,FileName{f}));
            writeDatFile(raw, FileName{f}, PathName);            
        end
        
    case 2
        % usage:  mat2dat(logsout, 'test_data.dat') 
        % params: 
        %           logsout:  Simulink.SimulationData.Dataset type data stored via simulation output
        %           'test_data.dat' :  filename to convert data to dat
        % returns: 
        %           'test_data.dat':  data converted to MDF format
        % simulation model call; assume logsout variable in the workspace
        
        dataStream = varargin{1};
        filename=varargin{2};
        if isa(dataStream, 'Simulink.SimulationData.Dataset')
            rawTSData = classSimDataset(dataStream);
            writeDatFile(rawTSData, filename, strcat(pwd,'\'));
        else
            disp('ERROR:  Pass Simulation dataset values only')
        end
end
end

function writeDatFile(rawStructure, FileName, PathName)
%% pass raw structure to preprocessing  
dataStruct = classStruct(rawStructure);
try
    dataTable = struct2table(dataStruct);
catch Error
    disp('    ')
    disp(' *** ERROR IN MAT FILE *** ');
    disp(Error.message);
    disp(' --> Check MAT file, all signals must have same length');
    return
end

%% Writing to a file 
if strcmpi(FileName(end-2:end), 'dat')
    fileoutname = FileName;
else
    fileoutname = strcat(FileName(1:end-4),'_expMDA.dat');
end
disp(['Writing data to dat file -- "' fileoutname '"']) ;
fileoutname = strcat(PathName,fileoutname);
expMDF(dataTable, fileoutname);
disp('   ');
disp(['Formatted data written to ',fileoutname]);
disp(repmat(char('-'), 1, 50));
disp('  ');
end

%% converts timeseries to structure
function dout = classStruct(din)
dout =struct;
names = fieldnames(din);
for n = 1:length(names)
    data_inside = eval(strcat('din.',names{n}));
    if isa(data_inside, 'timeseries')
        if ~isfield(dout,'time')
            dout.time = single(data_inside.Time);
        end
        varout = strcat('dout.',names{n});
        eval([ varout '= single(data_inside.Data);']);
        
    elseif isstruct(data_inside)
        d = classStruct(data_inside);
        name_d = fieldnames(d);
        for nd = 1:length(name_d)
            varout = strcat('dout.',name_d{nd});
            varval = eval(strcat('d.',name_d{nd}));
            eval([ varout '= varval;']);        
        end
        
    elseif isa(data_inside, 'Simulink.SimulationData.Dataset')
        dset = classSimDataset(data_inside);
        d = classStruct(dset);
        name_d = fieldnames(d);
        for nd = 1:length(name_d)
            varout = strcat('dout.',name_d{nd});
            varval = eval(strcat('d.',name_d{nd}));
            eval([ varout '= varval;']);        
        end
    
    elseif isfloat(data_inside)
        if ~isfield(dout,'time')
            dout.time = single([0:1:length(data_inside)-1]');
        end
        varout = strcat('dout.',names{n});
        eval([ varout '= single(data_inside);']);
    end
end
end

%% converts Simulation Dataset to structure-timeseries
function uTimeSeries = classSimDataset(dataStream)
uTimeSeries = struct;
nuTimeSeries = struct;
for var = 1:numElements(dataStream)
    dataSignal = dataStream.getElement(var);
    if (isUniform(dataSignal.Values.TimeInfo)) && ...
        (dataSignal.Values.TimeInfo.Length > 2)
        eval(strcat('uTimeSeries.',dataSignal.Name, '= dataSignal.Values;'));
        eval('uts= dataSignal.Values;');
    else
        eval(strcat('nuTimeSeries.',dataSignal.Name, '= dataSignal.Values;'));
    end
end
if isempty(uts)
    error('No Uniform TimeVector Found: Try logging atleast one signal without the enable/tigger/function subsystem');
end
fnames = fieldnames(nuTimeSeries); 
for var = 1:length(fnames)
    varb = eval(strcat('nuTimeSeries.',fnames{var}));
    if varb.TimeInfo.Length == 1  % constant ports 
        tempTS = timeseries;
        tempTS.Name = varb.Name;
        tempTS.Time = uts.Time;
        tempTS.Data =  ones(length(uts.Time),1) .* varb.Data;
    else    % all other data formats
        if isfloat(varb.Data)
            Val = interp1( varb.Time, varb.Data, uts.Time, 'linear');
        else
            Val = floor(interp1( varb.Time, single(varb.Data), uts.Time, 'linear'));
        end
        Val(isnan(Val))=0.0;
        tempTS = timeseries;
        tempTS.Name = varb.Name;
        tempTS.Time = uts.Time;
        tempTS.Data = Val;
    end
    eval(strcat('uTimeSeries.',tempTS.Name, '= tempTS;'));
end
end