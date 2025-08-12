%% Quick analysis during experiment 
clear all
close all

% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-500; %removes spikes below this threshold
artefact_high=500; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=8; %How long after the trigger do you want to analyse spikes for(ms)? 

%% Load Intan parameters
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

dName='amplifier';
vFID = fopen([filepath filesep dName '.dat'],'r');

mem_check = dir([filepath filesep 'amplifier.dat']);
T = 256;
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * 30000);

% Blank artifact
denoiseIntan_sab(filepath, dName, T, 0, Startpoint_analyse, Overall_time_to_analyse);
trig = loadTrig(0);

allExtract_sab_1(dName,filepath,T,0,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

% Spikes align
d = Depth_s(0);
load("Exp1_Sim_Xia.sp.mat")
nSamples = size(sp{find(~cellfun(@isempty,sp),1)},2)-1;
time_axis = (0:nSamples-1)*(1000/FS); % time axis in ms
pre_ms = 0;      % time before trigger (ms)
post_ms = 30;    % time after trigger (ms)

% plot
for i= 1:nChn
   if isempty(sp{i}), continue; end

    thisChn = d(i);
    sp_data = sp{thisChn};
    sp_times = sp_data(:,1);  % spike time in ms
    sp_waveforms = sp_data(:,2:end); % spike waveforms

    aligned_spikes = [];

    for t = 1:length(trig)
        t0 = trig(t) / (FS/1000); % trigger time in ms 
        idx = find(sp_times >= t0 - pre_ms & sp_times <= t0 + post_ms);
        if ~isempty(idx)
            aligned_spikes = [aligned_spikes; sp_waveforms(idx,:)]; 
        end
    end

    spikes{i} = aligned_spikes;
end

figure();
for ch = 1:nChn
    subplot(4, 8, ch);
    hold on;
    spikes_ch = spikes{ch};
    ymax = max(abs(spikes_ch),[],"all");
    for tr = 1:min(size(spikes_ch,1),200)
        plot(time_axis, spikes_ch(tr,:),'k');
    end
    hold off;
    title(sprintf('Ch %d', ch));
    xlabel('Time (ms)');
    ylabel('µV');
    ylim([-500,500]);

end
