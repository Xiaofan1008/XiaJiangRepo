fileDIR = dir('*_exp_datafile_*.mat');
if isempty(fileDIR)
    return
end
fileDIR = fileDIR.name;
load(fileDIR,'StimParams');
