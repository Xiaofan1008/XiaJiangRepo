clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================================================================
%%                     USER SETTINGS
%% ================================================================

data_folder      = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Single3';
raster_chn_start = 1;
raster_chn_end   = 32;
Electrode_Type   = 1;    % 0: rigid, 1: flex, 2: 4-shank flex

% ---------------- RASTER PARAMETERS ----------------
ras_win         = [-20 50];   % ms window for raster
bin_ms_raster   = 1;           % bin size (ms)
smooth_ms       = 3;           % PSTH smoothing window (Gaussian)


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
    v  = unique(idx_all(rr));
    v  = v(v>0)';
    stimChPerTrial_all{t} = v;
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb,~,combClass] = unique(comb,'rows');
combClass_win = combClass;
nSets = size(uniqueComb,1);

% Pulse period decode
pulseTrain_all = cell2mat(StimParams(2:end,9));
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);

[PulsePeriods,~,pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

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

    for si = 1:nSets
        stimIdx   = uniqueComb(si,:);
        stimIdx   = stimIdx(stimIdx>0);
        setLabel  = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');

        for pi = 1:n_PULSE
            pulse_val = PulsePeriods(pi);
            trials_this = find(pulseIdx == pi & combClass_win == si);
            if isempty(trials_this), continue; end

            figure('Color','w','Name',sprintf('Ch %d | %s', ich, setLabel));
            tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

            ax1 = nexttile([3 1]); hold(ax1,'on');
            ax2 = nexttile;        hold(ax2,'on');

            psth_curves = cell(1, n_AMP);
            maxRate = 0;
            y_cursor = 0;
            yticks_vals = {};
            yticks_labels = {};

            for ai = 1:n_AMP
                color   = cmap(ai,:);
                amp_val = Amps(ai);

                amp_trials = find(ampIdx == ai & pulseIdx == pi & combClass_win == si);
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
            xlim(ax1, ras_win);
            ylim(ax1, [0 y_cursor]);

            yt = cell2mat(yticks_vals);
            yticks(ax1, yt);
            yticklabels(ax1, yticks_labels);

            title(ax1, sprintf('Raster — Ch %d | Set %s', ich, setLabel));

            % ---------------- Finalize PSTH ----------------
            for ai = 1:n_AMP
                plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.5);
            end

            xline(ax2,0,'r','LineWidth',1.4);
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