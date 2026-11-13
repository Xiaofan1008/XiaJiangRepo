%% ============================================================
%   Recording Setup Figure Generator: Stacked Trace Viewer
%   - Plots raw/filtered traces in separate rows (Stacked)
%   - Mimics a "Computer Monitor" display of recording channels
%   - Randomly selects 1 trial per folder to show stability over time
% ============================================================
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
% [MODIFIED 1] Replaced a single string with a cell array of folders
% Add the exact paths for your Hour 0, Hour 2, Hour 4, etc. datasets here:
data_folders = {
    '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim5', ... % Session 1 (e.g., Hour 0)
    '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim2'  ... % Session 2 (Replace with Hour 2 path)
    '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim3'  ...
    '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim6' 
    % Add more folders as needed...
};

% -- Selection --
channels_to_plot = [27];                % Choose ONE channel to show multiple trials for
amps_to_plot     = [5];                 % Choose ONE amplitude (uA)
ptd_to_plot      = [0];                 % Choose ONE PTD (ms)
sets_to_plot     = [1];                 % Choose ONE Set

% -- Display --
% [MODIFIED 2] Time window updated to [-10 10] to perfectly match your sample image
plot_window_ms   = [-10 20];            % Time window (ms)
y_scale_uv       = [-300 300];          % Y-axis limits (uV) for each row
scale_bar_uv     = 100;                 % Size of the scale bar (uV)

Electrode_Type   = 1;                   
% 'raw'  = amplifier.dat (Raw wideband)
% 'dn'   = amplifier_dn_sab.dat (Artifact removed)
% 'mu'   = <base_name>.mu_sab.dat (Filtered 300-6000Hz)
trace_type = 'mu';    

%% ====================== INITIALIZE FIGURE ======================
nFolders = length(data_folders);

% [MODIFIED 3] Force the figure to be a perfect IEEE square (8.89 x 8.89 cm)
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto');
t = tiledlayout(nFolders, 1, 'TileSpacing', 'none', 'Padding', 'compact');

% Get depth map once
d = Depth_s(Electrode_Type);
ch_plot = channels_to_plot(1);
ch_intan = d(ch_plot);

%% ====================== MULTI-FOLDER PLOTTING LOOP ======================
% [MODIFIED 4] The entire loading and plotting process is now wrapped in a loop
for f_idx = 1:nFolders
    
    curr_folder = data_folders{f_idx};
    if ~isfolder(curr_folder)
        fprintf('Folder not found: %s\n', curr_folder);
        nexttile; axis off; title('Folder Not Found');
        continue;
    end
    cd(curr_folder);
    
    %% ====================== BASE NAME ======================
    parts = split(curr_folder, filesep);
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
    
    %% ====================== READ HEADER & OPEN DATA ======================
    [amp_channels, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
    nChn = numel(amp_channels);
    
    fid = fopen(data_file,'r');
    if fid < 0
        fprintf('Cannot open %s in folder %d\n', data_file, f_idx);
        nexttile; axis off;
        continue;
    end
    
    %% ====================== LOAD METADATA ======================
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
    trig = loadTrig(0);
    fDIR = dir('*_exp_datafile_*.mat');
    S = load(fDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
    StimParams = S.StimParams; simN = S.simultaneous_stim; E_MAP = S.E_MAP;
    
    trialAmps = cell2mat(StimParams(2:end,16)); trialAmps = trialAmps(1:simN:end);
    postTrig = cell2mat(StimParams(2:end,6)); postTrig = postTrig(2:simN:end);
    
    stimNames = StimParams(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP(2:end));
    comb = zeros(S.n_Trials, simN);
    for t = 1:S.n_Trials, v = idx_all((t-1)*simN + (1:simN)); v=v(v>0); comb(t,1:numel(v)) = v(:)'; end
    [~,~,combClass] = unique(comb,'rows','stable');
    
    %% ====================== TIME SAMPLES ======================
    samp_win = round(plot_window_ms/1000 * FS);
    Nsamp = samp_win(2)-samp_win(1)+1;
    time_ms = (samp_win(1):samp_win(2)) / FS * 1000;
    
    %% ====================== FIND RANDOM TRIAL ======================
    ptd_val = ptd_to_plot * 1000; 
    set_id  = sets_to_plot(1);
    amp_val = amps_to_plot(1);
    
    tidx_set = find(combClass == set_id);
    tidx_ptd = find(ismember(postTrig, ptd_val));
    tidx_amp = find(trialAmps == amp_val);
    tlist = intersect(intersect(tidx_set, tidx_ptd), tidx_amp);
    
    ax = nexttile;
    
    if isempty(tlist)
        fprintf('No valid trials found in folder %d.\n', f_idx);
        axis off; title('No Trials Found');
    else
        % [MODIFIED 5] Pick exactly ONE random trial from this folder
        n_available = numel(tlist);
        rand_idx = randperm(n_available, 1); 
        tr = tlist(rand_idx);
        
        % Read Data for that specific trial
        start_idx = trig(tr) + samp_win(1);
        byte_pos  = start_idx * nChn * 2;
        fseek(fid, byte_pos, 'bof');
        data_block = fread(fid, [nChn, Nsamp], 'int16') * 0.195;
        trace = data_block(ch_intan,:);
        
        % Plot the trace
        plot(time_ms, trace, 'k', 'LineWidth', 1); hold on;
        
        % [MODIFIED 6] Moved red trigger line to t=0 exactly
        xline(0, 'r-', 'LineWidth', 1); 
        
        % Formatting
        ylim(y_scale_uv);
        xlim(plot_window_ms);
        box off;
        set(gca, 'YColor', 'none'); 
        
        % [MODIFIED 7] Added Scale Bar logic ONLY for the bottom tile
        if f_idx < nFolders
            set(gca, 'XColor', 'none'); % Hide X-axis for all upper traces
        else
            % We are on the last row! Draw the scale bar and X-axis.
            xlabel('Time (ms)', 'FontSize', 9);
            set(gca, 'XColor', 'k', 'LineWidth', 1.2, 'FontSize', 9);
            
            % Scale bar math: draw a vertical line 1 ms from the left edge
            x_bar = plot_window_ms(1) + 1; 
            % Start the bar near the bottom limit (e.g., -200 uV)
            y_bar_bottom = y_scale_uv(1) + 50; 
            y_bar_top = y_bar_bottom + scale_bar_uv;
            
            % Plot the literal black line
            plot([x_bar, x_bar], [y_bar_bottom, y_bar_top], 'k-', 'LineWidth', 2.0);
            
            % Add the "100 uV" text just to the right of the line
            text(x_bar + 0.3, y_bar_bottom + (scale_bar_uv/2), sprintf('%d \\muV', scale_bar_uv), ...
                'VerticalAlignment', 'middle', 'FontSize', 9);
        end
    end
    fclose(fid);
end