data = dir('*_exp_datafile_*.mat');
if isempty(data)
    n_REP = 1;
    return
end
data = data.name;
load(data,'n_REP');
%load(data,'n_REP_true');
%n_REP=n_REP_true;