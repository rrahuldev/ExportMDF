function mat2dat(varargin)
% The script is invoked two ways 
% 1: with no arguments to the script, in which case, the user selects the
% one or more mat files to convert to dat file format.  The output is the
% same name as matfile but appended with "_expMDF.dat"
% 
% 2: inline with other scripts or model where the syntax below is followed. 
%       usage:  MAT2DAT(logsout, 'test_data.dat') 
%       params: 
%           logsout:  Simulink.SimulationData.Dataset type data stored via simulation output
%           'test_data.dat' :  filename to convert data to dat
%       returns: 
%           'test_data.dat':  data converted to MDF format
%       simulation model call; assume logsout variable in the workspace

% Author: Rahul Rajampeta
% Updated: April 11, 2019

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
    fileoutname = strcat(FileName(1:end-4),'_expMDF.dat');
end
disp(' ');disp(['Writing data to dat file -- "' fileoutname '"']) ;
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
function allSignals = classSimDataset(dataStream)
% todo check if data is fixed step or variable step
% todo check isnan or isinf
% for fixed or variable step, timestamp should be the same

% 1: extract all the signals from logsout
allSignals = struct;
for var = 1:numElements(dataStream)
    dataSignal = dataStream.getElement(var);
    % check if signal value type is timeseries (i.e. signals or mux) 
    if isobject(dataSignal.Values)
        allSignals = UpdateSignalsList(allSignals, dataSignal);
    
    elseif isstruct(dataSignal.Values)     % it must be a bus type signal
        fnames = fieldnames(dataSignal.Values);
        tempSignal = dataSignal;
        for i=1:length(fnames)
            eval(['tempSignal.Values = dataSignal.Values.' fnames{i} ';']);
            tempSignal.Name = [dataSignal.Name '_' tempSignal.Values.Name];
            tempSignal.Values.Name = tempSignal.Name;
            allSignals = UpdateSignalsList(allSignals, tempSignal);            
        end
            
    else % unknown type
        % todo
    end
        
end

% 2: get the common timestamp for all the signals
commonTime = [];
signalNames = fieldnames(allSignals);
for i=1:length(signalNames)
    tempTS = eval(['allSignals.' signalNames{i}]);
    if tempTS.TimeInfo.Length > length(commonTime)
        commonTime = tempTS.Time;
    end
end
% first element the common timestamp is more than 2s, then set warning
if isempty(commonTime(1) > 2) 
    warning('No Uniform TimeVector Found: Try logging atleast one signal without the enable/tigger/function subsystem');
end

% 3: fix the non-uniform signals
helptext = true;
promptinput = true;
for i=1:length(signalNames)
    tempTS = eval(['allSignals.' signalNames{i}]);
    if (tempTS.TimeInfo.Length ~= length(commonTime))
        % display one time help text
        if (helptext)
            helptext = false;
            fprintf('\n');
            disp('---------------------------------------------------------------------------');
            fprintf('Some signals are enable/function/trigger subsystem signals and have missing values,\nchoose one of the option below for respective signals - \n');
            fprintf('\t "z" - fill with zeros (default)\n')
            fprintf('\t "h" - hold last sample value \n')
            fprintf('\t "i" - interpolate in between\n')
            fprintf('\t "za" - fill with zeros and use this option for rest of the signals\n')
            fprintf('\t "ha" - hold last sample values and use this option for rest of the signals\n')
            fprintf('\t "ia" - interpolate in between and use this option for rest of the signals\n')
            disp('---------------------------------------------------------------------------');
            fprintf('\n');
        end
        
        % prompt for fill type, default zero
        if (promptinput)
            commandwindow;
            reply = input([tempTS.Name ' - z/h/i/za/ha/ia [z]:'],'s');
            % any other input default to zero
            if (isempty(reply) ||...
                isempty(cell2mat(strfind({'z','h','i','za','ha','ia'}, lower(reply)))))
                reply = 'z';
            end
            % use same option for rest of the signals
            if length(reply)>1
                promptinput = false;
            end
        end
        
        % signal processing
        [~, idx]  = intersect(commonTime, tempTS.Time, 'stable');
        newData = zeros(size(commonTime));
        if strcmpi(reply(1),'z')
            % fill missing values with zeros
            newData(idx) = tempTS.Data;
            if ~promptinput
                fprintf([tempTS.Name ': zero fill\n']);
            end
            
        elseif strcmpi(reply(1),'h')
            % hold last samples
            newData(idx) = tempTS.Data - [0;tempTS.Data(1:end-1)];
            newData = cumsum(newData);
            if ~promptinput
                fprintf([tempTS.Name ': hold last value\n']);
            end
            
        elseif strcmpi(reply(1),'i')
            % interpolate missing values
            if (tempTS.TimeInfo.Length > 2)
                newData = interp1(tempTS.Time, tempTS.Data, commonTime, 'linear','extrap');
                if ~promptinput
                    fprintf([tempTS.Name ': interpolate in between\n']);
                end
            else
                fprintf('\n');
                fprintf(['"' tempTS.Name '" - signal has less than 3 samples to interpolate, defaulting to zero fill\n']);
                newData(idx) = tempTS.Data;
            end

        end
        newTS = timeseries(newData, commonTime, 'name', tempTS.Name);
        eval(['allSignals.' tempTS.Name ' =  newTS;']);

    end
end

end

function allSignals = UpdateSignalsList(allSignals, varSignal)
% check vectors vs signal
if (size(varSignal.Values.Data, 2) == 1) % only signal
    % fix time to round off to 5 digits
    varSignal.Values.Time = round(varSignal.Values.Time*1e5)/1e5;
    eval(['allSignals.' varSignal.Name ' = varSignal.Values;']);
else % vector
    tempTS  = varSignal.Values;
    % fix time to round off to 5 digits
    tempTS.Time = round(tempTS.Time*1e5)/1e5;
    for i = 1:size(varSignal.Values.Data, 2)
        tempTS.Data = varSignal.Values.Data(:,i);
        tempTS.Name = [varSignal.Values.Name '_signal' num2str(i)];
        eval(['allSignals.' tempTS.Name ' = tempTS;']);
    end
end
end

