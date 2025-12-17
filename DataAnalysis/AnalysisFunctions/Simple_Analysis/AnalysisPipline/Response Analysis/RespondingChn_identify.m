%% ========================================================================
%  Responding Channel Detection (Refined, no R = Responding... reads)
% ========================================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER INPUT =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq2_5ms';

Electrode_Type = 1;    % 0 = rigid, 1 = single-shank flex, 2 = four-shank flex

% ---- Choose detection mode ----
Detection_Mode = 1;     % 1 = FR rule, 2 = statistical test

% ---- Window definitions ----
baseline_win_ms = [-90 -10];        % baseline always this window
post_win_ms     = [2 20];           % Mode 1 post-stim window
baseDur_s = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;
postDur_s = (post_win_ms(2)      - post_win_ms(1))      / 1000;

% Mode 2: per-pulse window
post_win_singlepulse_ms = 5;        % e.g. 5 ms after EACH pulse

% Simple rule parameters (Mode 1)
k_SD = 3;                            % ≥ μ + k·σ
% Extra robustness thresholds (for BOTH baseline=0 and noisy cases)
min_total_baseline_spikes   = 0;    % across all trials
min_total_post_spikes       = 5;    % across all trials
min_frac_trials_with_spikes = 0.2;  % ≥20% trials must have ≥1 spike
min_abs_post_FR             = 5;    % mean post FR must be ≥ 5 sp/s

% Mode 2 filters
min_FR_post = 2;                     % minimum post FR to call responsive
min_frac_trials = 0.2;               % ≥20% trials must contain spikes

FS = 30000;                          % sampling rate (Hz)

%% ================= PREP & LOAD DATA =================
cd(data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder,'_');
if numel(u)>=4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];

if isfile(ssd_file)
    S = load(ssd_file);
    if     isfield(S,'sp_corr'), sp = S.sp_corr;
    elseif isfield(S,'sp_SSD'),  sp = S.sp_SSD;
    elseif isfield(S,'sp_in'),   sp = S.sp_in;
    else, error("SSD file missing usable spike variable."); end

elseif isfile(base_file)
    S = load(base_file);
    if isfield(S,'sp_clipped'), sp = S.sp_clipped;
    else, sp = S.sp; end
else
    error("No spike file found.");
end
nCh = numel(sp);

%% ================= LOAD BAD CHANNELS & BAD TRIALS =================
bad_file = [base_name '.BadChannels.mat'];
if isfile(bad_file)
    BadCh_perSet = load(bad_file).BadCh_perSet;   % cell{si}
else
    BadCh_perSet = {};
end

bad_trials_file = [base_name '.BadTrials.mat'];
if isfile(bad_trials_file)
    BadTrials = load(bad_trials_file).BadTrials;  % BadTrials{ich}
else
    BadTrials = cell(nCh,1);
end

%% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig    = loadTrig(0);
trig_ms = trig / FS * 1000;

%% ================= LOAD StimParams =================
fileDIR = dir('*_exp_datafile_*.mat');
S = load(fileDIR(1).name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams = S.StimParams;
sim_stim   = S.simultaneous_stim;
E_MAP      = S.E_MAP;
n_Trials   = S.n_Trials;

% ---- Amplitudes ----
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:sim_stim:end);
[Amps,~,ampIdx] = unique(trialAmps);
Amps(Amps==-1) = 0;
nAMP = numel(Amps);

% ---- PTDs ----
if sim_stim > 1
    PTD_all = cell2mat(StimParams(3:sim_stim:end,6)); % µs
else
    PTD_all = zeros(n_Trials,1);                      % single pulse
end
[PTDs,~,ptdIdx] = unique(PTD_all);
nPTD = numel(PTDs);

% ---- Stim sets (order-sensitive) ----
stimNames = StimParams(2:end,1);
[~,idx_all] = ismember(stimNames, E_MAP(2:end));

stimSeq = zeros(n_Trials, sim_stim);
for t = 1:n_Trials
    rr = (t-1)*sim_stim + (1:sim_stim);
    v = idx_all(rr); 
    v = v(v>0);
    stimSeq(t,1:numel(v)) = v;
end
[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

% ---- Depth map ----
d = Depth_s(Electrode_Type);

%% ======================================================================
%                   MAIN LOOP: SET × AMP × PTD × CHANNEL
%% ======================================================================

Responding = struct();
fprintf("\nResponding Channel Detection...\n");

for si = 1:nSets

    if si <= numel(BadCh_perSet)
        BadCh_thisSet = BadCh_perSet{si};
    else
        BadCh_thisSet = [];
    end

    Responding.set(si).stimChannels = uniqueComb(si, uniqueComb(si,:)>0);

    for ai = 1:nAMP
    for pi = 1:nPTD

        % Store amplitude & PTD information (once per condition)
        Responding.set(si).amp(ai).amp_value              = Amps(ai);
        Responding.set(si).amp(ai).ptd(pi).PTD_us         = PTDs(pi);
        Responding.set(si).amp(ai).ptd(pi).PTD_ms         = PTDs(pi)/1000;
        Responding.set(si).amp(ai).ptd(pi).amp_value      = Amps(ai);   % redundant but handy
        Responding.set(si).amp(ai).ptd(pi).set_index      = si;
        Responding.set(si).amp(ai).ptd(pi).amp_index      = ai;
        Responding.set(si).amp(ai).ptd(pi).ptd_index      = pi;

        base_trials = find(combClass==si & ampIdx==ai & ptdIdx==pi);

        for ich = 1:length(d)
            ch = d(ich);

            % Skip bad channels or empty channels
            if ismember(ich, BadCh_thisSet) || isempty(sp{ch})
                Responding.set(si).amp(ai).ptd(pi).channel(ich).is_responsive = false;
                continue;
            end

            % Per-channel bad trials
            if ich <= numel(BadTrials)
                bad_tr_ch = BadTrials{ich};
            else
                bad_tr_ch = [];
            end

            trials_this = setdiff(base_trials, bad_tr_ch);
            if isempty(trials_this)
                Responding.set(si).amp(ai).ptd(pi).channel(ich).is_responsive = false;
                continue;
            end

            sp_times = sp{ch}(:,1);

            % ================= BASELINE =================
            nTr = numel(trials_this);
            FR_baseline = zeros(1,nTr);
            baseWinDur  = (baseline_win_ms(2)-baseline_win_ms(1))/1000;

            for k = 1:nTr
                tr = trials_this(k);
                t0 = trig_ms(tr);
                mask = sp_times >= (t0 + baseline_win_ms(1)) & ...
                                sp_times <  (t0 + baseline_win_ms(2));
                FR_baseline(k) = sum(mask) / baseWinDur;
            end

            % mu_b = mean(FR_baseline);
            % sd_b = std(FR_baseline);
            mu_b = median(FR_baseline);
            MAD = median(abs(FR_baseline - mu_b));   % robust SD (MAD)
            sd_b = MAD / 0.6745;
            % ================= MODE 1 — MOdified Z-score rule =================           
            if Detection_Mode == 1
                % --- post-stim FR per trial ---
                FR_post = zeros(1,numel(trials_this));
                for k = 1:numel(trials_this)
                    tr = trials_this(k);
                    t0 = trig_ms(tr);
                    mask = sp_times >= (t0 + post_win_ms(1)) & ...
                                    sp_times <  (t0 + post_win_ms(2));
                    FR_post(k) = sum(mask) / postDur_s;
                end
            
                % --- robust baseline stats (you already saw median helps) ---
                mu_b = median(FR_baseline);
                sd_b = median(abs(FR_baseline - mu_b));   % robust SD (MAD style)
            
                % --- mean post FR ---
                mu_p = mean(FR_post);
            
                % --- total spike counts (baseline & post) ---
                total_baseline_spikes = sum(FR_baseline) * baseDur_s;
                total_post_spikes     = sum(FR_post)     * postDur_s;
            
                % --- fraction of trials with at least one post spike ---
                frac_post_trials = mean(FR_post > 0);
            
                % ---- baseline+SD rule ----
                pass_FR_rule = (mu_p >= mu_b + k_SD * sd_b);
            
                % ---- additional gates for false positives ----
                pass_baseline_spikes = (total_baseline_spikes >= min_total_baseline_spikes);
                pass_post_spikes     = (total_post_spikes     >= min_total_post_spikes);
                pass_frac_trials     = (frac_post_trials      >= min_frac_trials_with_spikes);
                pass_abs_FR          = (mu_p                  >= min_abs_post_FR);
            
                isResp = pass_FR_rule & pass_baseline_spikes & ...
                         pass_post_spikes & pass_frac_trials & pass_abs_FR;
            
                % ---- store everything (so you can inspect later) ----
                Responding.set(si).amp(ai).ptd(pi).channel(ich).mean_baseline_FR   = mu_b;
                Responding.set(si).amp(ai).ptd(pi).channel(ich).sd_baseline_FR     = sd_b;
                Responding.set(si).amp(ai).ptd(pi).channel(ich).mean_post_FR       = mu_p;
                Responding.set(si).amp(ai).ptd(pi).channel(ich).FR_baseline_all    = FR_baseline;
                Responding.set(si).amp(ai).ptd(pi).channel(ich).FR_post_all        = FR_post;     
                Responding.set(si).amp(ai).ptd(pi).channel(ich).total_baseline_spikes = total_baseline_spikes;
                Responding.set(si).amp(ai).ptd(pi).channel(ich).total_post_spikes     = total_post_spikes;          
                Responding.set(si).amp(ai).ptd(pi).channel(ich).frac_post_trials   = frac_post_trials;
                Responding.set(si).amp(ai).ptd(pi).channel(ich).is_responsive      = isResp;          
                continue;
            end

            % ================= MODE 2 — Statistical =================
            FR_post = [];
            PTD_ms  = PTDs(pi)/1000;

            for k = 1:nTr
                tr = trials_this(k);
                t0 = trig_ms(tr);

                % Pulse 1
                mask = sp_times >= t0 & ...
                                sp_times < (t0+post_win_singlepulse_ms);
                FR_post(end+1) = sum(mask) / (post_win_singlepulse_ms/1000);

                % Pulse 2 (if sequential)
                if PTD_ms > 0
                    mask = sp_times >= (t0+PTD_ms) & ...
                                    sp_times < (t0+PTD_ms+post_win_singlepulse_ms);
                    FR_post(end+1) = sum(mask) / (post_win_singlepulse_ms/1000);
                end
            end

            p = ranksum(FR_post, FR_baseline, 'tail','right');
            frac_nonzero = mean(FR_post > 0);
            mu_p = mean(FR_post);

            isResp = (p < 0.01) && ...
                     (mu_p >= min_FR_post) && ...
                     (frac_nonzero >= min_frac_trials);

            % ---- DIRECT STRUCT ASSIGNMENT ----
            Responding.set(si).amp(ai).ptd(pi).channel(ich).is_responsive    = isResp;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).mean_baseline_FR = mu_b;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).sd_baseline_FR   = sd_b;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).mean_post_FR     = mu_p;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).FR_baseline_all  = FR_baseline;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).FR_post_all      = FR_post;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).p                = p;
            Responding.set(si).amp(ai).ptd(pi).channel(ich).frac_nonzero     = frac_nonzero;

        end % ich
    end % pi
    end % ai
end % si

%% ================= SAVE RESULT =================
outfile = sprintf('%s_RespondingChannels.mat', base_name);
save(outfile, 'Responding', ...
    'Detection_Mode','baseline_win_ms','post_win_ms','post_win_singlepulse_ms', ...
    'k_SD','min_FR_post','min_frac_trials','Amps','PTDs');

fprintf("\nSaved responding-channel results → %s\n", outfile);