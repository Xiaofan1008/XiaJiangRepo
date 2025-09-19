local = strsplit(filepath,'\');
local = local{end};
local = local(1:end-14);
if ~exist([savepath '\PSTH\' local  '\'],'dir')
    mkdir([savepath '\PSTH\' local  '\']);
end
saveas(PSTH,[savepath '\PSTH\' local  '\amp_' thisAmp '_dur_' thisDur '_chn_' thisChn '.png']);
pause(0.010);
close(PSTH);
if ~exist([savepath '\COLORPLOT\' local  '\'],'dir')
    mkdir([savepath '\COLORPLOT\' local  '\']);
end
saveas(COLORPLOT,[savepath '\COLORPLOT\' local  '\amp_' thisAmp '_dur_' thisDur '_chn_' thisChn '.png']);
pause(0.010);
close(COLORPLOT);