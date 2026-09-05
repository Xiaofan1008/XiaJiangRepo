data = dir('*_exp_datafile_*.mat');
if isempty(data)
    n_Trials = 1;
    return
end
data = data.name;
load(data,'n_Trials');
