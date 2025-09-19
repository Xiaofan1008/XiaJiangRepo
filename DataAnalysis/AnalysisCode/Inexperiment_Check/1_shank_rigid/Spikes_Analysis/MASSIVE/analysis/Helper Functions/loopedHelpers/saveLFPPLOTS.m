local = strsplit(filepath,'\');
local = local{end};
local = local(1:end-14);
if ~exist([savepath '\LFP\' local  '\'],'dir')
    mkdir([savepath '\LFP\' local  '\']);
end
saveas(LFP,[savepath '\LFP\' local  '\amp_' thisAmp '_dur_' thisDur '_chn_' thisChn '.png']);
pause(0.010);
close(LFP);