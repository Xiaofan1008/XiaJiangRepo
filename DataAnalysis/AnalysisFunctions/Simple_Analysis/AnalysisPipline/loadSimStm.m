data = dir('*_exp_datafile_*.mat');
if isempty(data)
    simultaneous_stim = 1;
    return
end
data = data.name;
load(data,'simultaneous_stim');