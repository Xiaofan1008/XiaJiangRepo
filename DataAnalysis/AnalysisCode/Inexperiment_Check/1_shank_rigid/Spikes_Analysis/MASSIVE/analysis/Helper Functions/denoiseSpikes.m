function sp = denoiseSpikes(sp,index)
if isempty(sp)
    return
end
chk = range(sp(:,2:end),2) > 500;
sp(chk,:) = [];
% Each recording is given a specific template and r2t (a 4*std of the zero
% condition template)
[template,r2t] = loadSpikeTemplate;
template = template{index}; r2t = r2t{index};
% Each spike gets a sum-of-squares
r2 = (sp(:,2:end) - template).^2;
% Each r2 much greater than the loaded r2t is rejected
chk = r2 > 16*r2t; chk = sum(chk,2);
sp(chk>0,:) = [];
end