local = strsplit(filepath,'\');
local = local{end};
local = local(1:end-14);
if ~exist([savepath '\SPIKES\' local  '\'],'dir')
    mkdir([savepath '\SPIKES\' local  '\']);
end
saveas(POST,[savepath '\SPIKES\' local  '\amp_' thisAmp '_dur_' thisDur '_chn_' thisChn '.png']);
pause(0.010);
close(POST);