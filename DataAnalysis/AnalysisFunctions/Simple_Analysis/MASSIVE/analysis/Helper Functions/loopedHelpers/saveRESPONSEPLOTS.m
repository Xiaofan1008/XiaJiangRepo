local = strsplit(filepath,'\');
local = local{end};
local = local(1:end-14);
if ~exist([savepath '\RESPONSE\' local  '\'],'dir')
    mkdir([savepath '\RESPONSE\' local  '\']);
end
saveas(stimResponse,[savepath '\RESPONSE\' local  '\chn_' thisChn '.png']);
pause(0.010);
close(stimResponse);