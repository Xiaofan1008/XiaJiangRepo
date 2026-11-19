%% Raw / Denoised / Filtered Trace Viewer
% For mixed-prefix stimulation files
%
% One figure = channel × set × amplitude × ISI
% Subplots = Prefix 1/2/3/4/5 + Simultaneous reference
%
% Main purpose:
%   Check whether the raw stimulation timing is correct.
%
% Important checks:
%   1) Prefix 1 should usually have only one active row.
%   2) Prefix 2 should usually have two active rows.
%   3) Prefix 5 should usually have five active rows.
%   4) If Prefix 1 already shows 5 artifacts, check PulseNum and PulsePeriod.

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1_260604_122628';

channels_to_plot = 1;        % channels to plot, using Depth_s index
amps_to_plot     = [10];     % amplitudes to include, [] means all
Electrode_Type   = 2;        % 0 rigid, 1 flex, 2 4-shank flex

% Prefix selection.
% [] means all detected prefixes.
prefix_to_plot  = [];

% ISI selection in ms.
% [] means all detected ISIs.
isi_to_plot_ms  = [];

% Set/order selection.
% [] means all detected condition set IDs.
sets_to_plot    = [];

% If 1, add full simultaneous condition in the same figure.
include_simultaneous_reference = 1;

% Number of raw traces to plot per condition.
nTrials_to_plot  = 30;

% Time window around trigger.
plot_window_ms   = [-5 25];

% Trace type:
%   'raw' = amplifier.dat
%   'dn'  = amplifier_dn_sab.dat
%   'mu'  = <base_name>.mu_sab.dat
trace_type = 'raw';

% Y-axis limit.
ylim_raw = [-7000, 7000];

% Print example StimParams rows for each prefix.
debug_print_trial_content = 1;

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder)
    error('Invalid folder');
end
cd(data_folder);

%% ====================== BASE NAME ======================
parts = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');

if numel(u) >= 4
    base_name = lastfld(1:u(end-1)-1);
else
    base_name = lastfld;
end

fprintf('\nData folder:\n%s\n', data_folder);
fprintf('Base name: %s\n', base_name);

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
        data_file  = [base_name '.mu_sab.dat'];

    otherwise
        error('Unknown trace_type "%s". Use ''raw'', ''dn'', or ''mu''.', trace_type);
end

%% ====================== READ INTAN HEADER ======================
fileinfo = dir('info.rhs');
[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);

fprintf('Sampling rate: %.1f Hz\n', FS);
fprintf('Number of amplifier channels: %d\n', nChn);

%% ====================== OPEN DATA FILE ======================
fid = fopen(data_file,'r');
if fid < 0
    error('Cannot open %s', data_file);
end

fprintf('Opened data file: %s\n', data_file);

%% ====================== LOAD TRIGGERS ======================
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD StimParams ======================
fDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fDIR),'No *_exp_datafile_*.mat found.');

S = load(fDIR(1).name);

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;
E_MAP             = S.E_MAP;

if isfield(S, 'CHN')
    CHN = S.CHN;
else
    CHN = [];
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match trigger number (%d). Using min of both for plotting.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== REMOVE HEADER ROW ======================
StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;
if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

%% ====================== TRIAL METADATA FROM StimParams COLUMNS ======================
% This is the important correction.
%
% Do not rely on separate arrays such as active_electrode_count_by_trial,
% because they may not be randomized together with StimParams.
%
% Instead, derive metadata directly from randomized StimParams columns:
%   26 = ActiveElectrodeCount
%   27 = PrefixLength
%   28 = ISI_ms
%   29 = ConditionType
%   30 = ConditionSetID

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30. Cannot use mixed-prefix metadata.');
end

firstRow_eachTrial = 1:simultaneous_stim:size(StimParams_data,1);

activeCount_trial    = cell2mat(StimParams_data(firstRow_eachTrial,26));
prefixLength_trial   = cell2mat(StimParams_data(firstRow_eachTrial,27));
isi_ms_trial         = cell2mat(StimParams_data(firstRow_eachTrial,28));
conditionType_trial  = cell2mat(StimParams_data(firstRow_eachTrial,29));
conditionSetID_trial = cell2mat(StimParams_data(firstRow_eachTrial,30));

fprintf('\nUsing trial metadata from StimParams columns 26–30.\n');

%% ====================== AMPLITUDES ======================
trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0 before unique().
trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);

if isempty(amps_to_plot)
    amps_sel = Amps;
else
    amps_sel = intersect(Amps, amps_to_plot);
end

%% ====================== LAST ACTIVE ARTIFACT TIME ======================
% Calculate final artifact time for each trial.
%
% For each active row:
%   final artifact time =
%       PTD_us + (PulseNum - 1) * PulsePeriod_us
%
% Then for the whole trial:
%   lastActivePTD_us(tr) = max(final artifact time across active rows)
%
% This also helps detect whether Prefix 1 has multiple pulses due to
% PulseNum > 1.

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials

    if conditionType_trial(tr) == 0
        lastActivePTD_us(tr) = 0;
        continue;
    end

    % For sequential prefix trials, prefix length is the number of active rows.
    if conditionType_trial(tr) == 1
        nActive_this = prefixLength_trial(tr);

    % For simultaneous trials, use active electrode count.
    elseif conditionType_trial(tr) == 2
        nActive_this = activeCount_trial(tr);

    else
        nActive_this = activeCount_trial(tr);
    end

    if isnan(nActive_this) || nActive_this < 1
        lastActivePTD_us(tr) = 0;
        continue;
    end

    nActive_this = min(round(nActive_this), simultaneous_stim);

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    activeRows = rr(1:nActive_this);

    % Column 6: post-trigger delay, us.
    ptd_us = cell2mat(StimParams_data(activeRows,6));

    % Column 8: pulse train number.
    pulseNum = cell2mat(StimParams_data(activeRows,8));

    % Column 9: pulse train period, us.
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    % Protect invalid values.
    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowLastArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowLastArtifact_us);
end

%% ====================== APPLY USER FILTERS ======================
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(sets_to_plot)
    set_sel = SetIDs_all;
else
    set_sel = intersect(SetIDs_all, sets_to_plot);
end

prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));

if isempty(prefix_to_plot)
    prefix_sel = prefix_all;
else
    prefix_sel = intersect(prefix_all, prefix_to_plot);
end

isi_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(isi_to_plot_ms)
    isi_sel = isi_all;
else
    isi_sel = intersect(isi_all, isi_to_plot_ms);
end

fprintf('\nDetected amplitudes: ');
disp(Amps');

fprintf('Selected amplitudes: ');
disp(amps_sel');

fprintf('\nDetected set IDs: ');
disp(SetIDs_all');

fprintf('Selected set IDs: ');
disp(set_sel');

fprintf('\nDetected prefixes: ');
disp(prefix_all');

fprintf('Selected prefixes: ');
disp(prefix_sel');

fprintf('\nDetected ISIs (ms): ');
disp(isi_all');

fprintf('Selected ISIs (ms): ');
disp(isi_sel');

fprintf('\nDetected condition types: ');
disp(unique(conditionType_trial)');

fprintf('Detected final active artifact times (us): ');
disp(unique(lastActivePTD_us)');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================
% This section prints one example trial for each prefix.
% It helps check whether the stimulation file is actually correct.

if debug_print_trial_content == 1 && ~isempty(set_sel) && ~isempty(amps_sel) && ~isempty(isi_sel)

    fprintf('\n================ DEBUG TRIAL CONTENT CHECK ================\n');

    debug_set = set_sel(1);
    debug_amp = amps_sel(1);
    debug_isi = isi_sel(1);

    fprintf('Debug Set = %d | Amp = %.1f uA | ISI = %.1f ms\n', ...
        debug_set, debug_amp, debug_isi);

    for ip = 1:numel(prefix_sel)

        prefix_val = prefix_sel(ip);

        tr_debug = find(conditionSetID_trial == debug_set & ...
                        conditionType_trial == 1 & ...
                        prefixLength_trial == prefix_val & ...
                        isi_ms_trial == debug_isi & ...
                        trialAmps == debug_amp, ...
                        1, 'first');

        if isempty(tr_debug)
            fprintf('\nPrefix %d: no matching trial found.\n', prefix_val);
            continue;
        end

        rr = (tr_debug-1)*simultaneous_stim + (1:simultaneous_stim);

        stimNames_debug = StimParams_data(rr,1);
        ptd_debug = cell2mat(StimParams_data(rr,6));
        pulseNum_debug = cell2mat(StimParams_data(rr,8));
        pulsePeriod_debug = cell2mat(StimParams_data(rr,9));
        amp_debug = cell2mat(StimParams_data(rr,16));
        activeCount_debug = cell2mat(StimParams_data(rr,26));
        prefix_debug = cell2mat(StimParams_data(rr,27));
        isi_debug = cell2mat(StimParams_data(rr,28));
        condType_debug = cell2mat(StimParams_data(rr,29));
        setID_debug = cell2mat(StimParams_data(rr,30));

        fprintf('\nPrefix %d | Trial %d\n', prefix_val, tr_debug);
        fprintf('  conditionType = %d, setID = %d, ISI = %.1f ms, amp = %.1f uA\n', ...
            conditionType_trial(tr_debug), conditionSetID_trial(tr_debug), ...
            isi_ms_trial(tr_debug), trialAmps(tr_debug));

        fprintf('  StimNames:       ');
        disp(stimNames_debug');

        fprintf('  PTD_us:          ');
        disp(ptd_debug');

        fprintf('  PulseNum:        ');
        disp(pulseNum_debug');

        fprintf('  PulsePeriod_us:  ');
        disp(pulsePeriod_debug');

        fprintf('  Amp_col16:       ');
        disp(amp_debug');

        fprintf('  ActiveCount_col26: ');
        disp(activeCount_debug');

        fprintf('  Prefix_col27:      ');
        disp(prefix_debug');

        fprintf('  ISI_col28:         ');
        disp(isi_debug');

        fprintf('  CondType_col29:    ');
        disp(condType_debug');

        fprintf('  SetID_col30:       ');
        disp(setID_debug');

        fprintf('  FinalArtifact_us = %.1f\n', lastActivePTD_us(tr_debug));
    end

    fprintf('===========================================================\n');
end

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
                subplot_labels = {};
                subplot_types  = {};
                subplot_values = [];

                % Prefix subplots.
                for ip = 1:length(prefix_sel)
                    subplot_labels{end+1} = sprintf('Prefix %d', prefix_sel(ip));
                    subplot_types{end+1}  = 'prefix';
                    subplot_values(end+1) = prefix_sel(ip);
                end

                % Simultaneous reference.
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
                figName = sprintf('%sTrace_Ch%d_Set%d_Amp%g_ISI%.1fms', ...
                                  data_label, ch_plot, set_id, amp_val, isi_val);

                figure('Name',figName,'Color','w','Position',[150 150 1500 850]);

                %% ----- Subplot loop -----
                for i_sub = 1:nSub

                    thisType  = subplot_types{i_sub};
                    thisValue = subplot_values(i_sub);

                    if strcmp(thisType, 'prefix')

                        prefix_val = thisValue;

                        % Sequential prefix condition.
                        tlist = find(conditionSetID_trial == set_id & ...
                                     conditionType_trial == 1 & ...
                                     prefixLength_trial == prefix_val & ...
                                     isi_ms_trial == isi_val & ...
                                     trialAmps == amp_val);

                        % Keep only trials that have triggers.
                        tlist = tlist(tlist <= nTrials_use);

                        lastPTD_this = unique(lastActivePTD_us(tlist));

                        if isempty(lastPTD_this)
                            lastPTD_text = 'NA';
                        else
                            lastPTD_text = sprintf('%.1f', max(lastPTD_this)/1000);
                        end

                        title_text = sprintf('Prefix %d | Final artifact %s ms', ...
                                             prefix_val, lastPTD_text);

                    elseif strcmp(thisType, 'sim')

                        % Full simultaneous condition.
                        % ISI is ignored for simultaneous.
                        tlist = find(conditionSetID_trial == set_id & ...
                                     conditionType_trial == 2 & ...
                                     trialAmps == amp_val);

                        tlist = tlist(tlist <= nTrials_use);

                        lastPTD_this = unique(lastActivePTD_us(tlist));

                        if isempty(lastPTD_this)
                            lastPTD_text = 'NA';
                        else
                            lastPTD_text = sprintf('%.1f', max(lastPTD_this)/1000);
                        end

                        title_text = sprintf('Simultaneous | Final artifact %s ms', lastPTD_text);

                    else
                        tlist = [];
                        title_text = 'Unknown';
                    end

                    subplot(nRows, nCols, i_sub);
                    hold on;

                    if isempty(tlist)
                        title(sprintf('%s\nNo trials', title_text), 'Interpreter','none');
                        xline(0,'r--');
                        xlabel('Time (ms)');
                        ylabel('uV');
                        ylim(ylim_raw);
                        xlim(plot_window_ms);
                        continue;
                    end

                    % Only draw first N trials to avoid overcrowding.
                    tlist = tlist(1:min(nTrials_to_plot, numel(tlist)));

                    %% ----- Plot trial traces -----
                    nPlotted = 0;

                    for k = 1:numel(tlist)

                        tr = tlist(k);

                        if tr > length(trig)
                            continue;
                        end

                        start_idx = trig(tr) + samp_win(1);

                        % Skip if the requested window starts before file beginning.
                        if start_idx < 0
                            continue;
                        end

                        byte_pos = start_idx * nChn * 2;

                        fseek(fid, byte_pos, 'bof');
                        data_block = fread(fid, [nChn, Nsamp], 'int16') * 0.195;

                        % Skip incomplete reads near the end of file.
                        if size(data_block,2) < Nsamp
                            continue;
                        end

                        plot(time_ms, data_block(ch_intan,:), 'LineWidth', 1);
                        nPlotted = nPlotted + 1;
                    end

                    %% ----- Axis formatting -----
                    xline(0,'r--');

                    % Draw final artifact line if it is after 0 ms.
                    if ~isempty(lastPTD_this)
                        finalArtifact_ms = max(lastPTD_this) / 1000;
                        if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                            xline(finalArtifact_ms, 'k--', 'LineWidth', 1);
                        end
                    end

                    xlabel('Time (ms)');
                    ylabel('uV');
                    ylim(ylim_raw);
                    xlim(plot_window_ms);

                    title(sprintf('%s (%d trials)', title_text, nPlotted), ...
                          'Interpreter','none');

                end % subplot loop

                %% ----- Figure title -----
                sgtitle(sprintf('%s Trace | Rec Ch %d | Set %d: %s | %g uA | Seq ISI %.1f ms', ...
                    data_label, ch_plot, set_id, set_label, amp_val, isi_val), ...
                    'Interpreter','none');

            end % ISI
        end % amplitude
    end % set
end % channel

fclose(fid);