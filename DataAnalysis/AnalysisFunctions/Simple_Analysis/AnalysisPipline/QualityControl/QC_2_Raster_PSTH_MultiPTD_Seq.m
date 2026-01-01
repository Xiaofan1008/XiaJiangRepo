clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================================================================
%%                     USER SETTINGS
%% ================================================================
data_folder      = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim1';
raster_chn_start = 1;
raster_chn_end   = 64;
Electrode_Type   = 2;    % 0: 1-rigid, 1: 1-flex, 2: 4-shank flex

% [NEW] Select specific PTDs to plot (e.g. [5 10]). Leave [] to plot ALL.
target_PTDs      = [5 ];   

% ---------------- RASTER PARAMETERS ----------------
ras_win         = [-50 80];   % ms window for raster
bin_ms_raster   = 1;           % bin size (ms)
smooth_ms       = 5;           % PSTH smoothing window (Gaussian)

% ================================================================
%                       PATH CHECK + BASE NAME
% ================================================================
if ~isfolder(data_folder)
    error('Specified folder does not exist.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

parts       = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder,'_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);
else
    base_name = last_folder;
end

% ================================================================
%                     LOAD FILTERED SPIKES
% ================================================================
fprintf('Loading filtered spikes from sp_xia_SSD.mat...\n');
ssd_file = [base_name '.sp_xia_SSD.mat'];
assert(isfile(ssd_file), ...
    'SSD Filtered spike file %s not found. Run QC_SpikeFilter.m first.', ssd_file);
S = load(ssd_file);
sp_clipped = S.sp_corr; 

% ================================================================
%                          LOAD TRIGGERS
% ================================================================
if isempty(dir('*.trig.dat'))
    cur_dir = pwd; 
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);
FS = 30000;

%% ================================================================
%%                     LOAD STIM PARAMS & DECODE
%% ================================================================
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No exp_datafile found.');
S = load(fileDIR(1).name, ...
    'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP             = S.E_MAP;
n_Trials          = S.n_Trials;

% Amplitude decode
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simultaneous_stim:end);
[Amps,~,ampIdx] = unique(trialAmps(:));
Amps(Amps==-1) = 0;
n_AMP = numel(Amps);
cmap  = lines(n_AMP); 

% Stim set decode
stimNames = StimParams(2:end,1);
E_NAME    = E_MAP(2:end);
[~, idx_all] = ismember(stimNames, E_NAME);

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    % KEEP ORDER (do NOT use unique)
    v = idx_all(rr);
    v = v(v>0);
    stimChPerTrial_all{t} = v(:)';     % row vector
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows', 'stable');
nSets = size(uniqueComb,1);
combClass_win = combClass;

% Pulse period decode
pulseTrain_all = cell2mat(StimParams(2:end,9));
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);
[PulsePeriods,~,pulseIdx] = unique(pulseTrain(:)); 
n_PULSE = numel(PulsePeriods);

% PTD Extraction (Pulse Train Delay)
% Column 6 usually contains PTD in microseconds. 
% For simultaneous_stim > 1, the delay is usually on the 2nd row of the pair.
PTD_all_us = cell2mat(StimParams(2:end, 6)); 
% We grab the PTD associated with the trial (usually the max delay in the group)
% Or specifically the delay of the second pulse.
trialPTD_us = PTD_all_us(2:simultaneous_stim:end); 
trialPTD_ms = trialPTD_us / 1000;

% Handle case where PTD might be 0 (Simultaneous)
% Create unique list of PTDs found in data
available_PTDs = unique(trialPTD_ms);

% Filter based on User Selection
if isempty(target_PTDs)
    use_PTDs = available_PTDs;
else
    % Keep only PTDs that exist in data AND are requested
    use_PTDs = intersect(available_PTDs, target_PTDs);
    if isempty(use_PTDs)
        warning('No requested PTDs found in dataset. Plotting nothing.');
    end
end
fprintf('Plotting PTDs: %s ms\n', num2str(use_PTDs'));

% Electrode map
d = Depth_s(Electrode_Type);

% ================================================================
%                   PREPARE RASTER/PSTH VARIABLES
% ================================================================
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;
g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

% ================================================================
%                         RASTER LOOP
% ================================================================
for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end
    S_ch = sp_clipped{ch};

    % Loop over Sets
    for si = 1:nSets
        stimIdx   = uniqueComb(si,:);
        stimIdx   = stimIdx(stimIdx>0);
        setLabel  = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');
        
        % Loop over Selected PTDs
        for pi = 1:length(use_PTDs)
            ptd_val = use_PTDs(pi);
            
            % Identify trials matching this Set AND this PTD
            % Note: We are ignoring PulsePeriods loop for simplicity unless needed, 
            % assuming PTD is the main variable here. If you have changing PulsePeriods 
            % AND PTDs, you need nested loops. 
            % Here we filter by PTD tolerance (float check).
            
            trials_this = find(combClass_win == si & abs(trialPTD_ms - ptd_val) < 0.001);
            
            if isempty(trials_this), continue; end
            
            % Define Title based on PTD
            if ptd_val == 0
                mode_str = 'Simultaneous';
            else
                mode_str = sprintf('Sequential (%g ms)', ptd_val);
            end
            
            figure('Color','w','Name',sprintf('Ch %d | %s | %s', ich, setLabel, mode_str));
            tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

            % Top 3 tiles for Raster, Bottom 1 for PSTH
            ax1 = nexttile([3 1]); hold(ax1,'on');
            ax2 = nexttile;        hold(ax2,'on');
            
            psth_curves = cell(1, n_AMP);
            maxRate = 0;
            y_cursor = 0;
            yticks_vals = {};
            yticks_labels = {};
            
            % Loop Amplitudes
            for ai = 1:n_AMP
                color   = cmap(ai,:);
                amp_val = Amps(ai);
                
                % Intersect trials_this (Set & PTD) with Amp
                amp_trials = trials_this(trialAmps(trials_this) == amp_val);
                
                nTr = numel(amp_trials);
                counts = zeros(1, numel(edges)-1);
                spkTimesPerTrial = cell(nTr,1);
                
                for t = 1:nTr
                    tr = amp_trials(t);
                    t0 = trig(tr)/FS*1000;
                    tt = S_ch(:,1);
                    tt = tt(tt >= t0+ras_win(1) & tt <= t0+ras_win(2)) - t0;
                    spkTimesPerTrial{t} = tt;
                    counts = counts + histcounts(tt, edges);
                end
                
                % ---------------- Raster ----------------
                for t = 1:nTr
                    tt = spkTimesPerTrial{t};
                    y0 = y_cursor + t;
                    for k = 1:numel(tt)
                        plot(ax1, [tt(k) tt(k)], [y0-0.4 y0+0.4], ...
                            'Color', color, 'LineWidth', 1.2);
                    end
                end
                
                if nTr > 0
                    yticks_vals{end+1}   = y_cursor + nTr/2;
                    yticks_labels{end+1} = sprintf('%d µA', amp_val);
                end
                y_cursor = y_cursor + max(nTr,1);
                
                % ---------------- PSTH ----------------
                if nTr > 0
                    rate = counts / (nTr * bin_s);
                    rate_s = filter(g,1,rate);
                    psth_curves{ai} = rate_s;
                    maxRate = max(maxRate, max(rate_s));
                else
                    psth_curves{ai} = zeros(1, numel(ctrs));
                end
            end
            
            % ---------------- Finalize Raster ----------------
            xline(ax1,0,'r','LineWidth',1.4);
            if ptd_val > 0
                xline(ax1, ptd_val, 'k--', 'LineWidth', 1.0); % Mark 2nd pulse
            end
            xlim(ax1, ras_win);
            ylim(ax1, [0 y_cursor]);
            yt = cell2mat(yticks_vals);
            yticks(ax1, yt);
            yticklabels(ax1, yticks_labels);
            title(ax1, sprintf('Raster — Ch %d | %s | %s', ich, setLabel, mode_str));
            
            % ---------------- Finalize PSTH ----------------
            for ai = 1:n_AMP
                plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.5);
            end
            xline(ax2,0,'r','LineWidth',1.4);
            if ptd_val > 0
                xline(ax2, ptd_val, 'k--', 'LineWidth', 1.0);
            end
            xlim(ax2, ras_win);
            if maxRate > 0
                ylim(ax2,[0 ceil(maxRate*1.1/10)*10]);
            else
                ylim(ax2,[0 1]);
            end
            xlabel(ax2,'Time (ms)');
            ylabel(ax2,'Rate (sp/s)');
            legend(ax2, arrayfun(@(a)sprintf('%.0f µA', a), Amps, 'UniformOutput', false), ...
                   'Box','off');
        end
    end
end