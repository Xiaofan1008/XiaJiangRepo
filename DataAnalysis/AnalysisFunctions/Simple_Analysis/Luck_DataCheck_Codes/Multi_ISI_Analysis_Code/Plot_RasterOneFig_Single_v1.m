%% =============================================================
% SINGLE-STIMULATION RASTER + PSTH FOR DATA CLEANING
%
% PURPOSE
%   Plot raster and PSTH results for a single-electrode stimulation dataset.
%
% IMPORTANT
%   - Uses sp_corr whenever available.
%   - Does not remove bad trials or bad channels.
%   - Can ignore an existing RespondingChannels file.
%   - Prints raster-row to global-trial-ID mappings.
%   - Maps stimulation hardware names to raster/depth channel numbers.
%   - Highlights the stimulation channel with a red background.
%   - Does not save or modify experimental data.
% =============================================================

clear;
close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ========================= USER SETTINGS ==============================

data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';

% Probe type:
%   0 = rigid single shank
%   1 = flexible single shank
%   2 = flexible four shank
Electrode_Type = 2;

% Depth/raster channels to display
raster_chn_start = 1;
raster_chn_end   = 64;

% Raster and PSTH window in ms
ras_win = [-50 80];

% PSTH settings
bin_ms_raster = 1;
smooth_ms     = 5;

% Empty means all available amplitudes
Plot_Amps = [];

% Single stimulation has PTD = 0
Plot_PTDs = 0;

% Ignore old responding-channel labels during initial cleaning
Show_Responding_Highlight = true;

% Print raster-row to global-trial-ID tables
Print_Trial_Mapping = true;

fig_position = [50 50 1600 900];

FS = 30000;

%% ======================== INITIAL SETUP ===============================

if ~isfolder(data_folder)
    error('Folder not found:\n%s',data_folder);
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder));

cd(data_folder);

fprintf('\n');
fprintf('============================================================\n');
fprintf('SINGLE-STIMULATION RASTER + PSTH\n');
fprintf('============================================================\n');
fprintf('Dataset folder:\n%s\n\n',data_folder);

%% ========================= FIND SPIKE FILE ============================

ssd_files  = dir('*.sp_xia_SSD.mat');
base_files = dir('*.sp_xia.mat');

if ~isempty(ssd_files)

    if numel(ssd_files) > 1
        warning('Multiple SSD spike files found. Using: %s', ...
            ssd_files(1).name);
    end

    spike_file = ssd_files(1).name;
    base_name  = erase(spike_file,'.sp_xia_SSD.mat');

elseif ~isempty(base_files)

    if numel(base_files) > 1
        warning('Multiple base spike files found. Using: %s', ...
            base_files(1).name);
    end

    spike_file = base_files(1).name;
    base_name  = erase(spike_file,'.sp_xia.mat');

else
    error('No *.sp_xia_SSD.mat or *.sp_xia.mat file was found.');
end

%% ========================== LOAD SPIKES ===============================

available_variables = who('-file',spike_file);

if ismember('sp_corr',available_variables)

    spike_variable = 'sp_corr';

elseif ismember('sp_pca',available_variables)

    spike_variable = 'sp_pca';
    warning('sp_corr is missing. Using sp_pca instead.');

elseif ismember('sp_SSD',available_variables)

    spike_variable = 'sp_SSD';
    warning('sp_corr is missing. Using sp_SSD instead.');

elseif ismember('sp_clipped',available_variables)

    spike_variable = 'sp_clipped';
    warning('sp_corr is missing. Using sp_clipped instead.');

elseif ismember('sp',available_variables)

    spike_variable = 'sp';
    warning('sp_corr is missing. Using sp instead.');

else
    error('No usable spike variable was found in:\n%s',spike_file);
end

SpikeLoad = load(spike_file,spike_variable);
sp = SpikeLoad.(spike_variable);

nCh = numel(sp);

fprintf('Spike file: %s\n',spike_file);
fprintf('Spike variable: %s\n',spike_variable);
fprintf('Recording channels: %d\n',nCh);

%% ==================== LOAD BAD CHANNEL INFORMATION ====================

BadCh_perSet = {};

badch_file = [base_name '.BadChannels.mat'];

if isfile(badch_file)

    BadChLoad = load(badch_file);

    if isfield(BadChLoad,'BadCh_perSet')
        BadCh_perSet = BadChLoad.BadCh_perSet;
        fprintf('Bad-channel file loaded: %s\n',badch_file);
    else
        warning(['Bad-channel file exists, but BadCh_perSet is missing. ' ...
                 'No channels will be marked as bad.']);
    end

else
    fprintf('No bad-channel file found. This is allowed during cleaning.\n');
end

%% ===================== LOAD BAD TRIAL INFORMATION =====================

BadTrialsPerCh = {};

badtr_file = [base_name '.BadTrials.mat'];

if isfile(badtr_file)

    BadTrialLoad = load(badtr_file);

    if isfield(BadTrialLoad,'BadTrials')
        BadTrialsPerCh = BadTrialLoad.BadTrials;
        fprintf('Bad-trial file loaded: %s\n',badtr_file);
    else
        warning(['Bad-trial file exists, but BadTrials is missing. ' ...
                 'No trials will be marked as bad.']);
    end

else
    fprintf('No bad-trial file found. All trials will be displayed.\n');
end

%% ==================== LOAD RESPONDING CHANNELS ========================

Resp = [];
hasResp = false;

resp_file = [base_name '_RespondingChannels.mat'];

if Show_Responding_Highlight

    if isfile(resp_file)

        RespLoad = load(resp_file);

        if isfield(RespLoad,'Responding')
            Resp = RespLoad.Responding;
            hasResp = true;
            fprintf('Responding-channel file loaded: %s\n',resp_file);
        else
            warning(['Responding-channel file exists, but Responding ' ...
                     'is missing.']);
        end

    else
        fprintf('No responding-channel file found.\n');
    end

else
    fprintf(['Responding-channel highlighting is disabled. ' ...
             'Existing labels will be ignored.\n']);
end

%% ========================== LOAD TRIGGERS =============================

trigger_files = dir('*.trig.dat');

if isempty(trigger_files)

    fprintf('No trigger file found. Running trigger cleaning...\n');

    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);

    trigger_files = dir('*.trig.dat');
end

if isempty(trigger_files)
    error('No *.trig.dat file could be found or created.');
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Triggers loaded: %d\n',nTrig);

%% ================= LOAD EXPERIMENT PARAMETERS =========================

experiment_files = dir('*_exp_datafile_*.mat');

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    warning('Multiple experiment files found. Using: %s', ...
        experiment_files(1).name);
end

ExpLoad = load(experiment_files(1).name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');

required_fields = {
    'StimParams'
    'simultaneous_stim'
    'E_MAP'
    'n_Trials'
};

for k = 1:numel(required_fields)

    if ~isfield(ExpLoad,required_fields{k})
        error('Experiment file is missing variable: %s', ...
            required_fields{k});
    end
end

StimParams = ExpLoad.StimParams;
sim_stim   = ExpLoad.simultaneous_stim;
E_MAP      = ExpLoad.E_MAP;
n_Trials   = ExpLoad.n_Trials;

fprintf('Experiment file: %s\n',experiment_files(1).name);
fprintf('Experiment trials: %d\n',n_Trials);
fprintf('Stimulation events per trial: %d\n',sim_stim);

if sim_stim ~= 1
    warning(['This dataset has simultaneous_stim = %d, not 1. ' ...
             'It may not be a single-stimulation dataset.'],sim_stim);
end

if nTrig < n_Trials

    error(['Only %d triggers were loaded for %d experiment trials. ' ...
           'The trigger file should be checked before plotting.'], ...
           nTrig,n_Trials);

elseif nTrig > n_Trials

    warning(['There are %d triggers but %d experiment trials. ' ...
             'Only experiment trial IDs will be used.'], ...
             nTrig,n_Trials);
end

%% ========================= DECODE AMPLITUDES ==========================

trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:sim_stim:end);

trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps);

if isempty(Plot_Amps)
    Plot_Amps = Amps(:).';
end

fprintf('Available amplitudes: %s uA\n',num2str(Amps(:).'));

%% ============================ DECODE PTD ===============================

if sim_stim > 1

    PTD_all = cell2mat(StimParams(3:sim_stim:end,6));

else

    PTD_all = zeros(n_Trials,1);
end

[PTDs,~,ptdIdx] = unique(PTD_all);
nPTD = numel(PTDs);

fprintf('Available PTDs: %s ms\n',num2str(PTDs(:).'/1000));

%% ==================== DECODE STIMULATION SETS =========================

stimNames = StimParams(2:end,1);

[~,idx_all] = ismember(stimNames,E_MAP(2:end));

stimSeq = zeros(n_Trials,sim_stim);

for trial_id = 1:n_Trials

    rows_this_trial = ...
        (trial_id-1)*sim_stim + (1:sim_stim);

    stimulation_indices = idx_all(rows_this_trial);
    stimulation_indices = ...
        stimulation_indices(stimulation_indices > 0);

    stimSeq(trial_id,1:numel(stimulation_indices)) = ...
        stimulation_indices;
end

[uniqueComb,~,combClass] = ...
    unique(stimSeq,'rows','stable');

nSets = size(uniqueComb,1);

fprintf('Detected stimulation sets: %d\n',nSets);

%% ========================== ELECTRODE MAP =============================

% Depth_s reads info.rhs from the current dataset folder
d = Depth_s(Electrode_Type);

if isempty(d)
    error('Depth_s returned an empty channel map.');
end

depth_range = ...
    raster_chn_start:min(raster_chn_end,numel(d));

if isempty(depth_range)
    error('The selected raster channel range is empty.');
end

nChPlot = numel(depth_range);

fprintf('Depth channels plotted: %d to %d\n', ...
    depth_range(1),depth_range(end));

%% ============= MAP STIMULATION TO RASTER CHANNELS ====================

% For each stimulation set, convert its hardware electrode names into the
% depth-channel numbers displayed in the raster.
StimDepthChannels = cell(nSets,1);
StimHardwareNames = cell(nSets,1);

fprintf('\nStimulation-channel mapping:\n');

for si = 1:nSets

    stim_map_indices = ...
        uniqueComb(si,uniqueComb(si,:) > 0);

    depth_channels_this_set = [];
    hardware_names_this_set = {};

    for k = 1:numel(stim_map_indices)

        map_position = stim_map_indices(k)+1;

        if map_position > numel(E_MAP)
            warning('Set %d contains an invalid E_MAP index.',si);
            continue;
        end

        hardware_name = ...
            char(string(E_MAP{map_position}));

        hardware_names_this_set{end+1} = ...
            hardware_name; %#ok<SAGROW>

        hardware_channel = ...
            electrode_name_to_spike_channel(hardware_name);

        % d maps:
        %   raster/depth channel -> spike/hardware channel
        depth_channel = find(d == hardware_channel);

        if isempty(depth_channel)

            warning(['Stimulation electrode %s could not be mapped ' ...
                     'to a raster/depth channel.'],hardware_name);

        elseif numel(depth_channel) > 1

            warning(['Stimulation electrode %s mapped to multiple ' ...
                     'depth channels. Using the first one.'], ...
                     hardware_name);

            depth_channels_this_set(end+1) = ...
                depth_channel(1); %#ok<SAGROW>

        else

            depth_channels_this_set(end+1) = ...
                depth_channel; %#ok<SAGROW>
        end

        if ~isempty(depth_channel)
            fprintf('  Set %d: %s -> Ch %d\n', ...
                si,hardware_name,depth_channel(1));
        end
    end

    StimDepthChannels{si} = ...
        unique(depth_channels_this_set,'stable');

    StimHardwareNames{si} = ...
        hardware_names_this_set;
end

%% =========================== PSTH SETUP ===============================

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs = edges(1:end-1)+diff(edges)/2;
bin_s = bin_ms_raster/1000;

% Preserve the smoothing convention from the original script
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g/sum(g);

%% ======================= MAIN CONDITION LOOPS =========================

for si = 1:nSets

    % These are the channel numbers displayed in the raster
    stimDepthThisSet = StimDepthChannels{si};

    if isempty(stimDepthThisSet)

        setLabel = 'Stim channel unmapped';

    else

        setLabel = strjoin( ...
            arrayfun(@(channel_number) ...
            sprintf('Ch%d',channel_number), ...
            stimDepthThisSet, ...
            'UniformOutput',false), ...
            ' + ');
    end

    % Bad depth channels associated with this stimulation set
    if ~isempty(BadCh_perSet) && si <= numel(BadCh_perSet)

        badDepthThisSet = BadCh_perSet{si};

        if isempty(badDepthThisSet)
            badDepthThisSet = [];
        end

    else

        badDepthThisSet = [];
    end

    for aVal = Plot_Amps

        ai = find(abs(Amps-aVal) < 1e-6,1);

        if isempty(ai)
            fprintf('Set %d | %g uA is unavailable. Skipping.\n', ...
                si,aVal);
            continue;
        end

        for pi = 1:nPTD

            PTD_us = PTDs(pi);
            PTD_ms = PTD_us/1000;

            if ~isempty(Plot_PTDs) && ...
                    ~any(abs(Plot_PTDs-PTD_ms) < 1e-6)
                continue;
            end

            trials_this = find( ...
                combClass == si & ...
                ampIdx == ai & ...
                ptdIdx == pi);

            if isempty(trials_this)
                continue;
            end

            %% -------- PRINT GLOBAL TRIAL-ID MAPPING ------------------

            if Print_Trial_Mapping

                fprintf('\n');
                fprintf('Set %d (%s) | Amp %g uA | PTD %g ms\n', ...
                    si,setLabel,aVal,PTD_ms);

                fprintf('Raster row -> global trial ID:\n');

                trial_mapping = table( ...
                    (1:numel(trials_this))', ...
                    trials_this(:), ...
                    'VariableNames', ...
                    {'RasterRow','GlobalTrialID'});

                disp(trial_mapping);
            end

            %% ---------------- STIMULATION MODE ------------------------

            if sim_stim == 1
                stimModeStr = 'Single';
            elseif abs(PTD_ms) < 1e-6
                stimModeStr = 'Simultaneous';
            else
                stimModeStr = 'Sequential';
            end

            %% ---------------- CREATE FIGURE ---------------------------

            figTitle = sprintf( ...
                'Set %d (%s) | Amp %.1f uA | PTD %.1f ms | %s', ...
                si,setLabel,aVal,PTD_ms,stimModeStr);

            figure( ...
                'Color','w', ...
                'Name',figTitle, ...
                'NumberTitle','off', ...
                'Position',fig_position);

            tiledlayout( ...
                'flow', ...
                'TileSpacing','compact', ...
                'Padding','compact');

            sgtitle(figTitle, ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'Interpreter','none');

            %% ---------------- CHANNEL LOOP ----------------------------

            for idxDepth = 1:nChPlot

                ich = depth_range(idxDepth);
                ch  = d(ich);

                % Is this depth/raster channel stimulated?
                isStimCh = ismember(ich,stimDepthThisSet);

                % Invalid map entries cannot be plotted
                if ch < 1 || ch > nCh
                    nexttile;
                    axis off;
                    continue;
                end

                % Preserve the subplot even if this channel has no spikes,
                % so that a silent stimulation channel can still be red.
                if isempty(sp{ch})
                    sp_times = [];
                else
                    sp_times = double(sp{ch}(:,1));
                end

                % Previously marked bad trials, for information only
                if ~isempty(BadTrialsPerCh) && ...
                        ich <= numel(BadTrialsPerCh)

                    badTr_ch = BadTrialsPerCh{ich};

                    if isempty(badTr_ch)
                        badTr_ch = [];
                    end

                else

                    badTr_ch = [];
                end

                %% ------------- COLLECT TRIAL SPIKES ------------------

                nTr = numel(trials_this);

                allTrialSpikes = cell(nTr,1);
                counts = zeros(1,numel(edges)-1);

                for ti = 1:nTr

                    trial_id = trials_this(ti);
                    trigger_ms = trig(trial_id)/FS*1000;

                    relative_spikes = sp_times( ...
                        sp_times >= trigger_ms+ras_win(1) & ...
                        sp_times <= trigger_ms+ras_win(2)) ...
                        - trigger_ms;

                    allTrialSpikes{ti} = relative_spikes;

                    counts = counts + ...
                        histcounts(relative_spikes,edges);
                end

                %% ---------------- CALCULATE PSTH ---------------------

                if ~any(counts)

                    rate_s = zeros(size(ctrs));

                else

                    rate = counts/(nTr*bin_s);
                    rate_s = filter(g,1,rate);
                end

                maxRate = max(rate_s);
                yMaxPSTH = max(50,ceil(maxRate*1.1/10)*10);

                %% ---------------- CREATE SUBPLOT ---------------------

                ax = nexttile;
                hold(ax,'on');

                %% -------- RESPONDING/BAD-CHANNEL STATUS -------------

                isResp = false;
                respTag = '';

                if hasResp && ...
                        si <= numel(Resp.set) && ...
                        ai <= numel(Resp.set(si).amp) && ...
                        pi <= numel(Resp.set(si).amp(ai).ptd) && ...
                        ich <= numel( ...
                            Resp.set(si).amp(ai).ptd(pi).channel)

                    R = Resp.set(si).amp(ai).ptd(pi).channel(ich);

                    if isfield(R,'is_responsive') && ...
                            R.is_responsive

                        isResp = true;
                        respTag = ' RESP';
                    end
                end

                isBadCh = ismember(ich,badDepthThisSet);

                %% ---------------- BACKGROUND COLORS ------------------

                % Color priority:
                %   stimulation + bad = orange
                %   stimulation       = strong pale red
                %   responding + bad  = pale orange
                %   responding        = pale pink
                %   bad               = gray

                if isStimCh && isBadCh

                    set(ax,'Color',[1.00 0.72 0.50]);

                elseif isStimCh

                    set(ax,'Color',[1.00 0.72 0.72]);

                elseif isResp && isBadCh

                    set(ax,'Color',[1.00 0.90 0.75]);

                elseif isResp

                    set(ax,'Color',[1.00 0.88 0.88]);

                elseif isBadCh

                    set(ax,'Color',[0.95 0.95 0.95]);
                end

                if isStimCh

                    ax.XColor = [0.70 0 0];
                    ax.YColor = [0.70 0 0];
                    ax.LineWidth = 1.8;

                elseif isResp

                    ax.XColor = [0.60 0 0];
                    ax.YColor = [0.60 0 0];
                    ax.LineWidth = 1.4;
                end

                %% ---------------- LEFT AXIS: PSTH --------------------

                yyaxis(ax,'left');

                if any(rate_s)
                    plot(ax,ctrs,rate_s,'LineWidth',1.4);
                end

                xlim(ax,ras_win);
                ylim(ax,[0 yMaxPSTH]);
                ylabel(ax,'Rate (sp/s)');

                %% ---------------- RIGHT AXIS: RASTER -----------------

                yyaxis(ax,'right');

                for ti = 1:nTr

                    relative_spikes = allTrialSpikes{ti};

                    if isempty(relative_spikes)
                        continue;
                    end

                    plot(ax, ...
                        relative_spikes, ...
                        ti*ones(size(relative_spikes)), ...
                        '.', ...
                        'Color',[0 0 0], ...
                        'MarkerSize',4);
                end

                ylim(ax,[0 nTr+1]);

                % Show a small number of raster-row labels
                nRasterTicks = min(nTr,6);

                if nRasterTicks > 0

                    rasterTicks = unique(round( ...
                        linspace(1,nTr,nRasterTicks)));

                    set(ax,'YTick',rasterTicks);
                end

                xline(ax,0,'r--','LineWidth',1);
                xlim(ax,ras_win);

                %% ---------------- TITLE AND LABELS -------------------

                chLabel = sprintf('Ch %d',ich);

                if isStimCh
                    chLabel = [chLabel ' STIM'];
                end

                if isBadCh
                    chLabel = [chLabel ' BAD'];
                end

                chLabel = [chLabel respTag];

                title_handle = title(ax, ...
                    chLabel, ...
                    'FontSize',11, ...
                    'FontWeight','bold', ...
                    'Interpreter','none');

                if isStimCh

                    set(title_handle, ...
                        'Color',[0.75 0 0], ...
                        'FontWeight','bold');

                elseif isResp

                    set(title_handle,'Color',[0.70 0 0]);
                end

                if idxDepth > ...
                        nChPlot-ceil(sqrt(nChPlot))
                    xlabel(ax,'Time (ms)');
                end

                %% -------- DISPLAY EXISTING BAD-TRIAL COUNT ----------

                if ~isempty(badTr_ch)

                    bad_count = numel( ...
                        intersect(trials_this,badTr_ch));

                    text(ax, ...
                        ras_win(1)+1, ...
                        nTr, ...
                        sprintf('badTr: %d',bad_count), ...
                        'FontSize',7, ...
                        'Color',[0.3 0.3 0.3], ...
                        'HorizontalAlignment','left', ...
                        'VerticalAlignment','top');
                end
            end
        end
    end
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('RASTER/PSTH PLOTTING COMPLETE\n');
fprintf('No trials or channels were removed.\n');
fprintf('No experimental files were modified.\n');
fprintf('============================================================\n');

%% ========================== LOCAL FUNCTION ============================

function spike_channel = electrode_name_to_spike_channel(electrode_name)
% Convert an electrode hardware name to the one-based spike-channel index
% used in sp_corr and Depth_s.
%
% Examples:
%   A-000 -> 1
%   A-021 -> 22
%   B-000 -> 33
%   C-000 -> 65
%   D-000 -> 97

electrode_name = upper(strtrim(char(electrode_name)));

tokens = regexp( ...
    electrode_name, ...
    '^([A-D])-(\d+)$', ...
    'tokens', ...
    'once');

if isempty(tokens)
    error('Invalid electrode name: %s',electrode_name);
end

bank_letter = tokens{1};
channel_zero_based = str2double(tokens{2});

switch bank_letter

    case 'A'
        bank_offset = 0;

    case 'B'
        bank_offset = 32;

    case 'C'
        bank_offset = 64;

    case 'D'
        bank_offset = 96;

    otherwise
        error('Unsupported electrode bank: %s',bank_letter);
end

spike_channel = ...
    bank_offset+channel_zero_based+1;

end