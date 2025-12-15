%% ============================================================
% SpikeFiltering_SSD_Corr.m
% Post-hoc spike waveform cleaning (SSD + template correlation)
% Input:  <base_name>.sp_xia.mat  (sp_clipped or sp)
% Output: <base_name>.sp_xia_SSD.mat (sp_in, sp_SSD, sp_corr)
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim6_251125_181554';

FS = 30000;                 % sampling frequency
SSD_threshold_factor = 16;  % SSD threshold (Allison-Walker style)
template_window_ms   = [2 50];   % window after trigger to build template
baseline_window_ms   = [-60 -5];  % baseline window for correlation filtering
corr_thresh          = 0.7;      % keep spikes with corr >= this in baseline

% filter choice
do_SSD_filter        = 1;
do_corr_filter       = 1;

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder)
    error('Folder does not exist.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

parts       = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1)-1);
else
    base_name = last_folder;
end

%% ================= LOAD SPIKES =================
fname_sp = [base_name '.sp_xia.mat'];
assert(isfile(fname_sp), 'Cannot find %s. Run your first spike-saving step first.', fname_sp);

S_in = load(fname_sp);
if isfield(S_in,'sp_clipped')
    sp_in = S_in.sp_clipped;
elseif isfield(S_in,'sp')
    sp_in = S_in.sp;
else
    error('No sp or sp_clipped found in %s', fname_sp);
end

nCh = numel(sp_in);

%% ================= LOAD TRIGGERS (for correlation filter) =================
if do_corr_filter
    if isempty(dir('*.trig.dat'))
        cur_dir = pwd; cleanTrig_sabquick; cd(cur_dir);
    end
    trig = loadTrig(0);
end

%% ================= 1) SSD-BASED FILTERING =================
sp_SSD = sp_in;

if do_SSD_filter
    fprintf('\nSSD-based waveform consistency filtering \n');

    for ch = 1:nCh
        if isempty(sp_in{ch}), continue; end

        waveforms = sp_in{ch}(:,2:end);  % exclude timestamp
        nSpikes   = size(waveforms,1);
        if nSpikes < 5
            sp_SSD{ch} = sp_in{ch};
            continue;
        end

        mean_wave = mean(waveforms,1);
        diff_mat  = waveforms - mean_wave;
        SSD       = sum(diff_mat.^2,2);
        mean_SSD  = mean(SSD);

        valid_idx = SSD <= (SSD_threshold_factor * mean_SSD);
        sp_SSD{ch} = sp_in{ch}(valid_idx,:);

        fprintf('Ch %2d: SSD kept %4d / %4d spikes (%.1f%%)\n', ...
            ch, sum(valid_idx), nSpikes, 100*sum(valid_idx)/nSpikes);
    end

    fprintf('SSD filtering complete.\n\n');
else
    fprintf('Skipping SSD filtering; sp_SSD = sp_in.\n');
end

%% ================= 2) TEMPLATE CORRELATION FILTERING =================
sp_corr = sp_SSD;

if do_corr_filter
    fprintf('\nTemplate-based correlation filtering\n');
    nTrials = numel(trig);

    for ch = 1:nCh
        if isempty(sp_corr{ch}), continue; end

        spt = sp_corr{ch}(:,1);     % times
        wfs = sp_corr{ch}(:,2:end); % waveforms
        nSpikes = size(wfs,1);
        if nSpikes < 10
            continue;
        end

        % ---- Build template from evoked spikes in template_window_ms ----
        evoked_mask = false(nSpikes,1);
        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            evoked_mask = evoked_mask | ...
                (spt >= t0 + template_window_ms(1) & ...
                 spt <= t0 + template_window_ms(2));
        end
        evoked_waves = wfs(evoked_mask,:);
        if size(evoked_waves,1) < 5
            fprintf('Ch %2d: not enough evoked spikes for template.\n', ch);
            continue;
        end
        template = mean(evoked_waves,1);

        % ---- Compute correlation of every spike vs template ----
        corr_vals = zeros(nSpikes,1);
        for i = 1:nSpikes
            corr_vals(i) = corr(template(:), wfs(i,:)');
        end

        % ---- Baseline trial-by-trial filtering ----
        final_keep = true(nSpikes,1);

        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            base_mask = (spt >= t0 + baseline_window_ms(1)) & ...
                        (spt <  t0 + baseline_window_ms(2));
            idx_trial = find(base_mask);

            if isempty(idx_trial), continue; end
            bad_idx = idx_trial(corr_vals(idx_trial) < corr_thresh);
            final_keep(bad_idx) = false;
        end

        kept = sum(final_keep);
        fprintf('Ch %2d: Corr filter kept %4d / %4d spikes (%.1f%%)\n', ...
            ch, kept, nSpikes, 100*kept/nSpikes);

        sp_corr{ch} = sp_corr{ch}(final_keep,:);
    end

    fprintf('Correlation filtering complete.\n');
else
    fprintf('Skipping correlation filtering; sp_corr = sp_SSD.\n');
end

%% ================= SAVE OUTPUT =================
QC_params = struct();
QC_params.SSD_threshold_factor = SSD_threshold_factor;
QC_params.template_window_ms   = template_window_ms;
QC_params.baseline_window_ms   = baseline_window_ms;
QC_params.corr_thresh          = corr_thresh;
QC_params.do_SSD_filter        = do_SSD_filter;
QC_params.do_corr_filter       = do_corr_filter;

out_name = [base_name '.sp_xia_SSD.mat'];
save(out_name, 'sp_in','sp_SSD','sp_corr','QC_params','-v7.3');
fprintf('\nSaved filtered spikes to %s\n', out_name);