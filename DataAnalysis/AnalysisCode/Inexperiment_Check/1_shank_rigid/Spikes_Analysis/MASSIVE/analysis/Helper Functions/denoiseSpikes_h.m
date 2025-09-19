function sp = denoiseSpikes_h(sp)
if isempty(sp)
    return
end
chk = range(sp(:,2:end),2) > 400;
sp(chk,:) = [];
chk = range(sp(:,2:end),2) < 115;
sp(chk,:) = [];
chk = diff(sp(:,1));
sp(chk < 1,:) = [];
chk = (sp(:,2) > -40 & sp(:,2) < 100);
sp(~chk,:) = [];
chk = (sp(:,9) > -60);
sp(~chk,:) = [];
chk = (sp(:,20) > -60);
sp(~chk,:) = [];
end