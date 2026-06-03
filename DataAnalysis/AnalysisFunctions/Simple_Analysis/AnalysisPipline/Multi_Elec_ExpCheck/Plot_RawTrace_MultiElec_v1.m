%% Raw / Denoised / Filtered Trace Viewer
% For mixed-prefix stimulation files
% One figure = channel × set × amplitude × ISI
% Subplots = Prefix 1/2/3/4/5 + Simultaneous reference

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

channels_to_plot = 15:20;                % channels to plot, using Depth_s index
amps_to_plot     = [10];                 % amplitudes to include, [] means all
Electrode_Type   = 2;                   % 0 rigid, 1 flex, 2 4-shank flex
% ================= XIA MODIFICATION =================
% For new mixed-prefix stimulation files:
%   Prefix 1 = first active electrode only
%   Prefix 2 = first two active electrodes
%   ...
%   Prefix 5 = full sequential 5-electrode condition
%
% isi_to_plot_ms is the actual ISI between sequential pulses.
% Example:
%   prefix 5 with ISI = 5 ms has artifacts at 0, 5, 10, 15, 20 ms.
prefix_to_plot  = [1 2 3 4 5];           % prefixes to plot, [] means all
isi_to_plot_ms  = [5];                   % ISI values in ms, [] means all
sets_to_plot    = [];                    % condition set IDs, [] means all

% If this is 1, the full simultaneous condition will be plotted together
% with the prefix subplots in the same figure.
include_simultaneous_reference = 1;
% =====================================================

nTrials_to_plot  = 30;                  % how many trials to draw per condition
plot_window_ms   = [-5 25];             % window around trigger
% 'raw'  = amplifier.dat
% 'dn'   = amplifier_dn_sab.dat   (denoised)
% 'mu'   = <base_name>.mu_sab.dat (filtered / MUA)
trace_type = 'raw';    % raw trace
% trace_type = 'dn';   % artifact blanked trace
% trace_type = 'mu';   % bandpass filtered trace

ylim_raw = [-7000, 7000];               % y-limit for raw traces

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder), error('Invalid folder'); end
cd(data_folder);

%% ====================== BASE NAME ======================
parts = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');
if numel(u)>=4
    base_name = lastfld(1:u(end-1)-1);      % e.g. 'Xia_Exp1_Single2'
else
    base_name = lastfld;
end

%% ====================== CHOOSE FILE BY TRACE TYPE ======================
switch trace_type
    case 'raw'
        data_label = 'Raw';
        data_file  = 'amplifier.dat';

    case 'dn'
        data_label = 'Denoised';
        data_file  = 'amplifier_dn_sab.dat';

    case 'mu'
        data_label = 'Filtered';
        data_file  = [base_name '.mu_sab.dat'];  % e.g. 'Xia_Exp1_Single2.mu_sab.dat'

    otherwise
        error('Unknown trace_type "%s". Use ''raw'', ''dn'', or ''mu''.', trace_type);
end

%% ====================== READ INTAN HEADER ======================
fileinfo = dir('info.rhs');
[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);

%% ====================== OPEN DATA FILE ======================
fid = fopen(data_file,'r');
if fid < 0
    error('Cannot open %s', data_file);
end

%% ====================== LOAD TRIGGERS ======================
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

%% ====================== LOAD StimParams ======================
fDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fDIR),'No *_exp_datafile_*.mat found.');

% ================= XIA MODIFICATION =================
% Load the full experiment file so we can use the new mixed-prefix metadata.
S = load(fDIR(1).name);

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;   % now means rows/slots per trial
n_Trials          = S.n_Trials;
E_MAP             = S.E_MAP;

if isfield(S, 'CHN')
    CHN = S.CHN;
else
    CHN = [];
end
% =====================================================

%% ====================== AMPLITUDES ======================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);

% Convert inactive/zero-control amplitude from -1 to 0 before unique().
trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);

if isempty(amps_to_plot)
    amps_sel = Amps;
else
    amps_sel = intersect(Amps, amps_to_plot);
end

%% ====================== LOAD NEW METADATA ======================
% ================= XIA MODIFICATION =================
% New mixed-prefix files should contain these variables.
% If they exist, use them directly.
%
% conditionType:
%   0 = zero-control
%   1 = sequential prefix/recovery
%   2 = full simultaneous
%
% conditionSetID:
%   which entered stimulation order this trial belongs to.
if isfield(S, 'active_electrode_count_by_trial')
    activeCount_trial    = S.active_electrode_count_by_trial(:);
    prefixLength_trial   = S.prefix_length_by_trial(:);
    isi_ms_trial         = S.isi_ms_by_trial(:);
    conditionType_trial  = S.condition_type_by_trial(:);
    conditionSetID_trial = S.condition_set_id_by_trial(:);
else
    % Backward-compatible fallback for old files.
    % Old files do not have prefix information, so treat every trial as one
    % stimulation condition with active count = simultaneous_stim.
    fprintf('\nNo mixed-prefix metadata found. Using old-style fallback.\n');

    activeCount_trial    = simultaneous_stim * ones(n_Trials,1);
    prefixLength_trial   = simultaneous_stim * ones(n_Trials,1);
    conditionType_trial  = ones(n_Trials,1);
    conditionSetID_trial = ones(n_Trials,1);

    % Estimate ISI/PTD from the last stored stimulation row.
    postTrig_all = cell2mat(StimParams(2:end,6));
    if simultaneous_stim > 4
        postTrig = postTrig_all(5:simultaneous_stim:end);
    elseif simultaneous_stim == 4
        postTrig = postTrig_all(4:simultaneous_stim:end);
    elseif simultaneous_stim == 3
        postTrig = postTrig_all(3:simultaneous_stim:end);
    elseif simultaneous_stim == 2
        postTrig = postTrig_all(2:simultaneous_stim:end);
    else
        postTrig = zeros(n_Trials,1);
    end
    isi_ms_trial = postTrig(:) ./ 1000;
end
% =====================================================

%% ====================== LAST ACTIVE PTD ======================
% ================= XIA MODIFICATION =================
% Calculate the final artifact timing for each trial.
% This is mainly used for title/checking.
% Example:
%   Prefix 3, ISI 5 ms -> lastActivePTD_us = 10000 us.
lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials
    activeCount_this = activeCount_trial(tr);

    if isnan(activeCount_this) || activeCount_this < 1
        lastActivePTD_us(tr) = 0;
    else
        activeCount_this = min(round(activeCount_this), simultaneous_stim);
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this; % +1 for StimParams header
        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end
% =====================================================

%% ====================== APPLY USER FILTERS ======================
% Select condition sets
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(sets_to_plot)
    set_sel = SetIDs_all;
else
    set_sel = intersect(SetIDs_all, sets_to_plot);
end

% Select prefixes
prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));

if isempty(prefix_to_plot)
    prefix_sel = prefix_all;
else
    prefix_sel = intersect(prefix_all, prefix_to_plot);
end

% Select ISIs for sequential prefix trials
isi_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(isi_to_plot_ms)
    isi_sel = isi_all;
else
    isi_sel = intersect(isi_all, isi_to_plot_ms);
end

fprintf('\nDetected amplitudes:'); disp(Amps');
fprintf('Selected amplitudes:'); disp(amps_sel');

fprintf('\nDetected set IDs:'); disp(SetIDs_all');
fprintf('Selected set IDs:'); disp(set_sel');

fprintf('\nDetected prefixes:'); disp(prefix_all');
fprintf('Selected prefixes:'); disp(prefix_sel');

fprintf('\nDetected ISIs (ms):'); disp(isi_all');
fprintf('Selected ISIs (ms):'); disp(isi_sel');

fprintf('\nDetected condition types:'); disp(unique(conditionType_trial)');
fprintf('Detected last active PTDs (us):'); disp(unique(lastActivePTD_us)');

%% ====================== TIME SAMPLES ======================
samp_win = round(plot_window_ms/1000 * FS);
Nsamp = samp_win(2)-samp_win(1)+1;
time_ms = (samp_win(1):samp_win(2)) / FS * 1000;

%% ====================== CHANNEL MAP ======================
d = Depth_s(Electrode_Type);

%% ====================== MAIN LOOP ======================
% Plot structure:
%   one figure = channel × set × amplitude × ISI
%   subplots = prefixes + optional simultaneous reference

for ich = 1:length(channels_to_plot)

    ch_plot  = channels_to_plot(ich);
    ch_intan = d(ch_plot);

    if ch_intan > nChn
        warning('Channel %d maps to Intan channel %d, larger than nChn=%d. Skipped.', ...
                ch_plot, ch_intan, nChn);
        continue;
    end

    for i_set = 1:length(set_sel)
        set_id = set_sel(i_set);

        %% ----- Build set label -----
        % For new files, conditionSetID corresponds to the row/order in CHN.
        if ~isempty(CHN) && set_id <= size(CHN,1)
            stimVec = CHN(set_id,:);
            stimVec = stimVec(stimVec > 0);
            set_label = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimVec,'UniformOutput',false),'→');
        else
            set_label = sprintf('Set%d', set_id);
        end

        for i_amp = 1:length(amps_sel)
            amp_val = amps_sel(i_amp);

            for i_isi = 1:length(isi_sel)
                isi_val = isi_sel(i_isi);

                %% ----- Decide subplot list -----
                % Prefix subplots
                subplot_labels = {};
                subplot_types  = {};   % 'prefix' or 'sim'
                subplot_values = [];   % prefix number, or 0 for simultaneous

                for ip = 1:length(prefix_sel)
                    subplot_labels{end+1} = sprintf('Prefix %d', prefix_sel(ip));
                    subplot_types{end+1}  = 'prefix';
                    subplot_values(end+1) = prefix_sel(ip);
                end

                % Add simultaneous reference into the same figure.
                if include_simultaneous_reference == 1
                    subplot_labels{end+1} = 'Simultaneous';
                    subplot_types{end+1}  = 'sim';
                    subplot_values(end+1) = 0;
                end

                nSub = length(subplot_labels);
                if nSub == 0
                    continue;
                end

                nRows = ceil(sqrt(nSub));
                nCols = ceil(nSub / nRows);

                %% ----- Create figure -----
                figName = sprintf('%sTrace_Ch%d_Set%d_Amp%d_ISI%.1fms', ...
                                  data_label, ch_plot, set_id, amp_val, isi_val);

                figure('Name',figName,'Color','w','Position',[150 150 1500 850]);

                %% ----- Subplot loop -----
                for i_sub = 1:nSub

                    thisType  = subplot_types{i_sub};
                    thisValue = subplot_values(i_sub);

                    if strcmp(thisType, 'prefix')
                        prefix_val = thisValue;

                        % Sequential prefix condition:
                        % same set, same amplitude, selected ISI, selected prefix.
                        tlist = find(conditionSetID_trial == set_id & ...
                                     conditionType_trial == 1 & ...
                                     prefixLength_trial == prefix_val & ...
                                     isi_ms_trial == isi_val & ...
                                     trialAmps == amp_val);

                        lastPTD_this = unique(lastActivePTD_us(tlist));
                        if isempty(lastPTD_this)
                            lastPTD_text = 'NA';
                        else
                            lastPTD_text = sprintf('%.1f', lastPTD_this(1)/1000);
                        end

                        title_text = sprintf('Prefix %d | Last PTD %s ms', ...
                                             prefix_val, lastPTD_text);

                    elseif strcmp(thisType, 'sim')
                        % Full simultaneous condition:
                        % include it as reference even though ISI_ms = 0.
                        tlist = find(conditionSetID_trial == set_id & ...
                                     conditionType_trial == 2 & ...
                                     trialAmps == amp_val);

                        title_text = sprintf('Simultaneous');

                    else
                        tlist = [];
                        title_text = 'Unknown';
                    end

                    subplot(nRows, nCols, i_sub); hold on

                    if isempty(tlist)
                        title(sprintf('%s\nNo trials', title_text), 'Interpreter','none');
                        xline(0,'r--');
                        xlabel('Time (ms)');
                        ylabel('µV');
                        ylim(ylim_raw);
                        continue;
                    end

                    % Only draw first N trials to avoid overcrowding.
                    tlist = tlist(1:min(nTrials_to_plot, numel(tlist)));

                    %% ----- Plot trial traces -----
                    for k = 1:numel(tlist)
                        tr = tlist(k);

                        if tr > length(trig)
                            continue;
                        end

                        start_idx = trig(tr) + samp_win(1);
                        byte_pos  = start_idx * nChn * 2;

                        % Skip if the requested window starts before file beginning.
                        if start_idx < 0
                            continue;
                        end

                        fseek(fid, byte_pos, 'bof');
                        data_block = fread(fid, [nChn, Nsamp], 'int16') * 0.195;

                        % Skip incomplete reads near the end of file.
                        if size(data_block,2) < Nsamp
                            continue;
                        end

                        plot(time_ms, data_block(ch_intan,:), 'LineWidth', 1);
                    end

                    %% ----- Axis formatting -----
                    xline(0,'r--');
                    xlabel('Time (ms)');
                    ylabel('µV');
                    ylim(ylim_raw);
                    xlim(plot_window_ms);

                    if strcmp(thisType, 'prefix')
                        title(sprintf('%s (%d trials)', title_text, numel(tlist)), ...
                              'Interpreter','none');
                    else
                        title(sprintf('%s (%d trials)', title_text, numel(tlist)), ...
                              'Interpreter','none');
                    end

                end % subplot loop

                %% ----- Figure title -----
                sgtitle(sprintf('%s Trace | Rec Ch %d | Set %d: %s | %d µA | Seq ISI %.1f ms', ...
                    data_label, ch_plot, set_id, set_label, amp_val, isi_val), ...
                    'Interpreter','none');

            end % ISI
        end % amplitude
    end % set
end % channel

fclose(fid);