function data = applyCAR(data)
if ~isempty(data)
    data = bsxfun(@minus, data, median(data,2)); % subtract median of each channel
    data = bsxfun(@minus, data, median(data,1)); % subtract median of each time point       
end