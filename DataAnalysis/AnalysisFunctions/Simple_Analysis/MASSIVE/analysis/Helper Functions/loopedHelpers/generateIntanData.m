%% Load in Intan data
if exist('info.rhs','file')
    filepath = pwd;
end
if ~exist('filepath','var')
    [~,filepath] = uigetfile('*.rhs','Select an Intan info file');
end
datafile = dir([filepath, '\*_exp_datafile*.mat']);
if isempty(datafile)
    disp('I need a datafile!');
    return;
end
datafile = ['\', datafile.name];
NN_cutoff = datetime('06-Nov-2018 00:00:00');
fileinfo = dir([filepath,'\info.rhs']);
if datetime(fileinfo.date) < NN_cutoff
    NN_order = 0;
else
    NN_order = 1;
end
if datetime(fileinfo.date) > datetime('28-Nov-2018 00:00:00')
    NN_I = 1;
else 
    NN_I = 0;
end
info = read_Intan_RHS2000_file(filepath,'\info.rhs',0);
FS = info.frequency_parameters.amplifier_sample_rate;       % Sampling frequency
nChn = length(info.amplifier_channels);

