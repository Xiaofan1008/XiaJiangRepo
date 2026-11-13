%% ============================================================
%   Recording Setup Figure Generator: Stacked Trace Viewer
%   - Plots raw/filtered traces in separate rows (Stacked)
%   - Mimics a "Computer Monitor" display of recording channels
%   - **Randomly selects** trials to plot each time
% ============================================================
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim3';

% -- Selection --
channels_to_plot = [27];                % Choose ONE channel to show multiple trials for
amps_to_plot     = [5];                % Choose ONE amplitude (uA)
ptd_to_plot      = [0];                 % Choose ONE PTD (ms)
sets_to_plot     = [1];                 % Choose ONE Set

% -- Display --
nTrials_to_plot  = 4;                   % Number of RANDOM trials to stack
plot_window_ms   = [-10 50];            % Time window (ms)
y_scale_uv       = [-300 300];          % Y-axis limits (uV) for each row

Electrode_Type   = 1;                   
% 'raw'  = amplifier.dat (Raw wideband)
% 'dn'   = amplifier_dn_sab.dat (Artifact removed)
% 'mu'   = <base_name>.mu_sab.dat (Filtered 300-6000Hz)
trace_type = 'mu';    

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder), error('Invalid folder'); end
cd(data_folder);

%% ====================== BASE NAME ======================
parts = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');
if numel(u)>=4, base_name = lastfld(1:u(end-1)-1); else, base_name = lastfld; end

%% ====================== CHOOSE FILE ======================
switch trace_type
    case 'raw', data_label = 'Raw'; data_file = 'amplifier.dat';
    case 'dn',  data_label = 'Denoised'; data_file = 'amplifier_dn_sab.dat';
    case 'mu',  data_label = 'Filtered'; data_file = [base_name '.mu_sab.dat'];
    otherwise, error('Unknown trace_type');
end

%% ====================== READ HEADER ======================
[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);

%% ====================== OPEN DATA ======================
fid = fopen(data_file,'r');
if fid < 0, error('Cannot open %s', data_file); end

%% ====================== LOAD METADATA ======================
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

fDIR = dir('*_exp_datafile_*.mat');
S = load(fDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams = S.StimParams; simN = S.simultaneous_stim; E_MAP = S.E_MAP;

% Amps
trialAmps = cell2mat(StimParams(2:end,16)); trialAmps = trialAmps(1:simN:end);
% PTDs
postTrig = cell2mat(StimParams(2:end,6)); postTrig = postTrig(2:simN:end);
% Sets
stimNames = StimParams(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP(2:end));
comb = zeros(S.n_Trials, simN);
for t = 1:S.n_Trials, v = idx_all((t-1)*simN + (1:simN)); v=v(v>0); comb(t,1:numel(v)) = v(:)'; end
[~,~,combClass] = unique(comb,'rows','stable');

%% ====================== TIME SAMPLES ======================
samp_win = round(plot_window_ms/1000 * FS);
Nsamp = samp_win(2)-samp_win(1)+1;
time_ms = (samp_win(1):samp_win(2)) / FS * 1000;
d = Depth_s(Electrode_Type);

%% ====================== PLOTTING LOOP ======================
% Filter user selection
ptd_val = ptd_to_plot * 1000; % Convert to us
set_id  = sets_to_plot(1);
amp_val = amps_to_plot(1);
ch_plot = channels_to_plot(1);
ch_intan = d(ch_plot);

% Find Trials
tidx_set = find(combClass == set_id);
tidx_ptd = find(ismember(postTrig, ptd_val));
tidx_amp = find(trialAmps == amp_val);
tlist = intersect(intersect(tidx_set, tidx_ptd), tidx_amp);

if isempty(tlist)
    fprintf('No trials found for this condition.\n');
else
    % --- RANDOM SELECTION ---
    n_available = numel(tlist);
    n_pick = min(nTrials_to_plot, n_available);
    
    if n_available > 0
        rand_indices = randperm(n_available, n_pick); % Pick random indices
        trials_to_show = tlist(sort(rand_indices));   % Sort to keep chronological order
    else
        trials_to_show = [];
    end
    
    % --- CREATE STACKED FIGURE ---
    figure('Color','w','Position',[200 200 600 150 * numel(trials_to_show)]);
    t = tiledlayout(numel(trials_to_show), 1, 'TileSpacing', 'none', 'Padding', 'compact');
    
    for k = 1:numel(trials_to_show)
        tr = trials_to_show(k);
        
        % Read Data
        start_idx = trig(tr) + samp_win(1);
        byte_pos  = start_idx * nChn * 2;
        fseek(fid, byte_pos, 'bof');
        data_block = fread(fid, [nChn, Nsamp], 'int16') * 0.195;
        trace = data_block(ch_intan,:);
        
        % Plot Tile
        ax = nexttile;
        plot(time_ms, trace, 'k', 'LineWidth', 1.2); hold on;
        xline(2.5, 'r-', 'LineWidth', 1.5); % Red trigger line
        
        % Formatting
        ylim(y_scale_uv);
        xlim(plot_window_ms);
        box off;
        set(gca, 'YColor', 'none'); % Hide Y axis line/ticks
        set(gca, 'XColor', 'none'); % This hides the bottom black line
        % Add Label on the right
        text(plot_window_ms(2), 0, sprintf('Trial %d', tr), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10, 'FontWeight', 'bold');
        
        % Only show X-axis for the bottom plot
        if k < numel(trials_to_show)
            set(gca, 'XColor', 'none');
        else
            xlabel('Time (ms)', 'FontWeight', 'bold');
            set(gca, 'XColor', 'k', 'LineWidth', 1.2);
        end
    end
    
    % sgtitle(sprintf('%s Trace | Ch %d | Set %d | %d µA', ...
    %     data_label, ch_plot, set_id, amp_val), 'FontWeight', 'bold');
end

fclose(fid);