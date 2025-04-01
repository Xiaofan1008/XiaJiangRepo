data = dir('*_exp_datafile_*.mat');
if isempty(data)
    return
end
data = data.name;
load(data,'StimParams');
