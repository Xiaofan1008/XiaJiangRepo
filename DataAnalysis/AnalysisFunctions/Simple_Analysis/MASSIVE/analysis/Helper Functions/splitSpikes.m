function splitSpikes
%% This function loads in an sp.mat file and splits it into stim and base sp files based on trig lines
disp('Splitting spikes');
sp = [];
name = dir('*.sp.mat');
name = name.name;
name = name(1:end-6);
loadSpikes;
nChn = size(sp,2);
loadTrig;
stim_sp = sp;
base_sp = sp;
for c = 1:nChn
    for t = 1:length(trig)
        chka = sp{c}(:,1) >= trig(t);
        chkb = sp{c}(:,1) <= trig(t) + 201;
        chk = and(chkb,chka);
        base_sp{c}(chk,:) = NaN;
    end
    stim_sp{c}(~isnan(base_sp{c}(:,1)),:) = [];
    base_sp{c}(isnan(base_sp{c}(:,1)),:) = [];
end
save([name 'basesp.mat'],'base_sp','-v7.3');
save([name 'stimsp.mat'],'stim_sp','-v7.3');
end