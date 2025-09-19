function threshold = loadThreshold
tmp = dir('threshold.mat');
if ~isempty(tmp)
    load(tmp.name, 'threshold')
else
    return;
end
end