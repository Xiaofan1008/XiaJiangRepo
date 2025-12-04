clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

pos_limit = 100;
neg_limit = -100;

%% User Params
spike_chn_start = 48;
spike_chn_end   = 48;
plot_amps       = [4 6];      % amplitudes to include
Electrode_Type  = 2;
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Single1_251203_121428';

%% Check folder
if ~isfolder(data_folder), error('Invalid folder'); end
cd(data_folder);

%% Extract base name
parts   = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');
if numel(u)>=4
    base_name = lastfld(1:u(end-1)-1);
else
    base_name = lastfld;
end

FS = 30000;

%% Load raw spike file
sp_files = dir(fullfile(data_folder,'*.sp.mat'));
assert(~isempty(sp_files),'No .sp.mat');
S = load(fullfile(data_folder, sp_files(1).name));
assert(isfield(S,'sp'),'No sp');
sp = S.sp;

%% Load triggers
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

%% Load clipped spike file (first pulse spikes)
% load([base_name '.sp_xia_FirstPulse.mat']);
% sp_clipped = sp_seq;

load([base_name '.sp_xia.mat']);

%% Load stim params
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR),'No exp file');
S = load(fileDIR(1).name,'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
n_Trials           = S.n_Trials;
E_MAP              = S.E_MAP;

%% Amplitudes
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simultaneous_stim:end);

[Amps,~,ampIdx] = unique(trialAmps(:));
Amps(Amps==-1) = 0;
n_AMP = numel(Amps);
cmap = lines(n_AMP);

ampMask = ismember(trialAmps, plot_amps);   % apply amplitude filtering

%% --------- Sequence-sensitive stimulation sets ---------
stimNames = StimParams(2:end,1);
[~,idx_all] = ismember(stimNames, E_MAP(2:end));

stimChPerTrial_all = cell(n_Trials,1);
for t=1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);

    % keep order, remove zeros only
    v = idx_all(rr);
    v = v(v > 0);

    stimChPerTrial_all{t} = v(:)';   % row vector
end

% build ordered combination matrix
comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t, 1:numel(v)) = v;
end

% FORCE sequence-sensitivity and stability
[uniqueComb,~,combClass] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);
combClass_win = combClass;

%% PTD decode
postTrigDelay_all = cell2mat(StimParams(2:end,6));
postTrigDelay = postTrigDelay_all(2:simultaneous_stim:end);

[uniqueDelays,~,delayIdx] = unique(postTrigDelay);
n_DELAYS = numel(uniqueDelays);

%% Preparations for plotting
d = Depth_s(Electrode_Type);
win_ms = 300;
bin_ms = 2;
nBins = 50/bin_ms;
amp_threshold = 100;
layout_row = 5;
layout_col = 5;

%% Main plotting
for ich = spike_chn_start:spike_chn_end

    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end

    sp_times = sp_clipped{ch}(:,1);
    sp_wave  = sp_clipped{ch}(:,2:end);

    ok = all(abs(sp_wave)<=amp_threshold,2);
    sp_times = sp_times(ok);
    sp_wave  = sp_wave(ok,:);
    if isempty(sp_times), continue; end

    t_wave = (0:size(sp_wave,2)-1)/FS*1000;

    for delay_i = 1:n_DELAYS
        delay_val = uniqueDelays(delay_i);

        for s = 1:nSets
            set_id = s;

            % SEQUENCE-SENSITIVE set + PTD + amplitude filter
            trial_mask = (combClass_win == set_id) & ...
                         (delayIdx == delay_i) & ...
                         ampMask;

            if ~any(trial_mask), continue; end

            trial_ids = find(trial_mask);

            stimIdx = uniqueComb(set_id,:);
            stimIdx = stimIdx(stimIdx>0);

            figName = sprintf('Record Ch%d | Stim Ch %s | Inter-Stimulus-Interval %d ms | Amps %s', ...
                    ich, strjoin(string(stimIdx),','), delay_val/1000, mat2str(plot_amps));
            figure('Name',figName,'Color','w','Position',[100 100 1400 800]);

            tiledlayout(layout_row,layout_col,'Padding','compact','TileSpacing','compact');

            all_spikes_by_bin_amp = cell(nBins,n_AMP);

            for idt = 1:numel(trial_ids)
                tr = trial_ids(idt);

                t0_ms = trig(tr)/FS*1000;
                amp_id = ampIdx(tr);

                mask_sp = sp_times>=t0_ms & sp_times<(t0_ms+win_ms);
                if ~any(mask_sp), continue; end

                rel_t = sp_times(mask_sp) - t0_ms;
                wavef = sp_wave(mask_sp,:);

                for j = 1:numel(rel_t)
                    b = floor(rel_t(j)/bin_ms)+1;
                    if b>=1 && b<=nBins
                        all_spikes_by_bin_amp{b,amp_id}(end+1,:) = wavef(j,:);
                    end
                end
            end

            allW = cell2mat(all_spikes_by_bin_amp(:));
            if isempty(allW), continue; end

            y_max = max(abs(allW(:)));
            y_lim = [-1 1]*ceil(y_max/50)*50;

            for b = 1:nBins
                nexttile; hold on

                for a = 1:n_AMP
                    if ~ismember(Amps(a), plot_amps), continue; end

                    waves = all_spikes_by_bin_amp{b,a};
                    if isempty(waves), continue; end

                    aligned = zeros(size(waves));
                    for k = 1:size(waves,1)
                        [~,midx] = min(waves(k,:));
                        shift = ceil(size(waves,2)/2)-midx;
                        aligned(k,:) = circshift(waves(k,:),shift,2);
                    end

                    plot(t_wave,aligned','Color',[cmap(a,:) 0.3]);
                end

                cnt = 0;
                for a = 1:n_AMP
                    if ismember(Amps(a), plot_amps)
                        cnt = cnt + size(all_spikes_by_bin_amp{b,a},1);
                    end
                end

                title(sprintf('%d-%d ms (%d spikes)', ...
                      (b-1)*bin_ms, b*bin_ms, cnt));
                xlabel('Time (ms)');
                ylabel('µV');
                ylim(y_lim);
                axis square;
            end

            % legend for chosen amplitudes
            keepA = find(ismember(Amps,plot_amps));
            h_leg = gobjects(numel(keepA),1);

            for ii = 1:numel(keepA)
                a = keepA(ii);
                h_leg(ii) = plot(nan,nan,'-','Color',cmap(a,:),'LineWidth',1.5);
            end

            legend(h_leg, arrayfun(@(x) sprintf('%g µA',x), ...
                Amps(keepA), 'UniformOutput',false), 'Location','northeastoutside');
        end
    end
end