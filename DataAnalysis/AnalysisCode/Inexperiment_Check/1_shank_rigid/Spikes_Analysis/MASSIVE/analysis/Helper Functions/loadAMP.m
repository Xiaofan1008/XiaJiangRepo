function AMP = loadAMP
tmp = dir('*exp_datafile_*.mat');
if ~isempty(tmp)
    load(tmp.name, 'AMP')   
else
    AMP = [];
    return;
end
end