data = dir('*_exp_datafile_*.mat');
if isempty(data)
    DUALSTIM = 1;
    return
end
data = data.name;
load(data,'DUALSTIM');
