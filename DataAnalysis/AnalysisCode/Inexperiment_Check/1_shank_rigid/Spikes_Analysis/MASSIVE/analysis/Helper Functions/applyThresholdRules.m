function [all,sig,nsig] = applyThresholdRules(tmp1,tmp2,tmp3,bR,r2lim,sat,maxfR,CJ)
all(:,1) = fieldnames(tmp1);
all(:,2) = struct2cell(tmp1);
sig(:,1) = fieldnames(tmp2);
sig(:,2) = struct2cell(tmp2);
nsig(:,1) = fieldnames(tmp3);
nsig(:,2) = struct2cell(tmp3);
% Rule 1: Baseline firing rate exceeds 5 spikes per second
ind = cell2mat(all(22,2)) < bR;
for i = 1:size(all,1)
    all{i,2}(ind) = [];
end
ind = cell2mat(sig(22,2)) < bR;
for i = 1:size(sig,1)
    sig{i,2}(ind) = [];
end
ind = cell2mat(nsig(22,2)) < bR;
for i = 1:size(nsig,1)
    nsig{i,2}(ind) = [];
end
% Rule 2: r^2 must exceed 0.9
ind = cell2mat(all(21,2)) < r2lim;
for i = 1:size(all,1)
    all{i,2}(ind) = [];
end
ind = cell2mat(sig(21,2)) < r2lim;
for i = 1:size(sig,1)
    sig{i,2}(ind) = [];
end
ind = cell2mat(nsig(21,2)) < r2lim;
for i = 1:size(nsig,1)
    nsig{i,2}(ind) = [];
end
% Rule 3: Channel must saturate
ind = cell2mat(all(25,2)) == sat;
for i = 1:size(all,1)
    all{i,2}(ind) = [];
end
ind = cell2mat(sig(25,2)) == sat;
for i = 1:size(sig,1)
    sig{i,2}(ind) = [];
end
ind = cell2mat(nsig(25,2)) == sat;
for i = 1:size(nsig,1)
    nsig{i,2}(ind) = [];
end
% Rule 4: Channel must reach 50 spikes per second
ind = cell2mat(all(23,2)) < maxfR;
for i = 1:size(all,1)
    all{i,2}(ind) = [];
end
ind = cell2mat(sig(23,2)) < maxfR;
for i = 1:size(sig,1)
    sig{i,2}(ind) = [];
end
ind = cell2mat(nsig(23,2)) < maxfR;
for i = 1:size(nsig,1)
    nsig{i,2}(ind) = [];
end
% Rule 5: Rodent or Marmoset Data?
ind = cell2mat(all(26,2)) ~= CJ;
for i = 1:size(all,1)
    all{i,2}(ind) = [];
end
ind = cell2mat(sig(26,2)) ~= CJ;
for i = 1:size(sig,1)
    sig{i,2}(ind) = [];
end
ind = cell2mat(nsig(26,2)) ~= CJ;
for i = 1:size(nsig,1)
    nsig{i,2}(ind) = [];
end
%% Convert back to structure
all = cell2struct(all(:,2),all(:,1));
sig = cell2struct(sig(:,2),sig(:,1));
nsig = cell2struct(nsig(:,2),nsig(:,1));
end