function DUR = loadDUR
tmp = dir('*exp_datafile_*.mat');
if ~isempty(tmp)
    load(tmp.name, 'DUR')
else
    return;
end
end