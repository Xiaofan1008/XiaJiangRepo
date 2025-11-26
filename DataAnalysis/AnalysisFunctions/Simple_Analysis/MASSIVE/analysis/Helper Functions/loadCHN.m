filepath = pwd;
tmp = dir([filepath filesep '*exp_datafile_*.mat']);
if ~isempty(tmp)
    tmp = [filepath filesep tmp.name];
    CHN = load(tmp,'CHN');
    CHN = CHN.CHN;
end