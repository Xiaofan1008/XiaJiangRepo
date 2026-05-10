%% ============================================================
%   Spatial Map of Responding Channels (Corrected Loading)
% ============================================================
%  - Loads 'RespondingChannels.mat' directly
%  - Plots electrode geometry (Red = Responsive, Gray = Non-Responsive)
%  - Handle Simultaneous OR Sequential folders automatically
% ============================================================
clear;

%% === User Inputs ===
% Point this to the folder you want to visualize
data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Single1_251125_110714'; 

amp_to_plot = 10;   % µA (Must match an amplitude used in the experiment)
Electrode_Type = 1; % 1: Single shank (Flex), 2: Four shank

% Geometry Settings
nShank = 1;
chPerShank = 32;
spacing_y = 1;   % Vertical spacing
spacing_x = 4;   % Distance between shanks (if multi-shank)

%% === 1. Load Experiment Data & Responding Logic ===
if ~isfolder(data_folder), error('Folder does not exist.'); end
cd(data_folder);

% A. Load Responding Logic
f_resp = dir('*RespondingChannels.mat');
if isempty(f_resp), error('No RespondingChannels.mat found!'); end
RespStruct = load(f_resp(1).name).Responding;
fprintf('Loaded: %s\n', f_resp(1).name);

% B. Load Experiment Parameters (for E_MAP and Amps)
f_exp = dir('*_exp_datafile_*.mat');
if isempty(f_exp), error('No Experiment Datafile found!'); end
S = load(f_exp(1).name, 'StimParams', 'E_MAP', 'simultaneous_stim','n_Trials');
Stim = S.StimParams;
E_MAP = S.E_MAP;
simN = S.simultaneous_stim;

% C. Detect Depth / Channel Mapping
d = Depth_s(Electrode_Type); 
nCh_total = length(d);

%% === 2. Parse Amplitudes & Sets ===
% We need to find which "Index" corresponds to 'amp_to_plot'
amps_all = cell2mat(Stim(2:end,16));
trialAmps = amps_all(1:simN:end);
[Amps_Unique, ~, ~] = unique(trialAmps);
Amps_Unique(Amps_Unique==-1) = 0;

amp_idx = find(abs(Amps_Unique - amp_to_plot) < 0.1, 1);
if isempty(amp_idx)
    error('Amplitude %.1f µA not found in this dataset. Available: %s', ...
        amp_to_plot, mat2str(Amps_Unique'));
end

fprintf('Target Amplitude: %.1f µA (Index %d)\n', amp_to_plot, amp_idx);

%% === 3. Assign Coordinates (Shank Geometry) ===
channel_groups = cell(1, nShank);
ch_counter = 1;
for s = 1:nShank
    start_idx = ch_counter;
    end_idx = min(start_idx + chPerShank - 1, nCh_total);
    channel_groups{s} = start_idx:end_idx;
    ch_counter = end_idx + 1;
end

x_pos = []; y_pos = [];
for s = 1:nShank
    nCh_shank = numel(channel_groups{s});
    y_vals = (1:nCh_shank) * spacing_y;
    x_vals = ones(1, nCh_shank) * (s-1) * spacing_x;
    x_pos = [x_pos, x_vals];
    y_pos = [y_pos, y_vals];
end

%% === 3b. Parse Stimulation Channels (New Addition) ===
% Extract stimulation names and map them to channel numbers
stimNames = Stim(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

% Create the combination matrix
nTr = S.n_Trials;
comb = zeros(nTr, simN);
for t = 1:nTr
    rr = (t-1)*simN + (1:simN);
    v = idx_all(rr); v = v(v>0);
    comb(t,1:numel(v)) = v(:).';
end

% Get unique sets (this aligns with sIdx 1, 2, etc.)
[uniqueComb, ~, ~] = unique(comb, 'rows', 'stable');

%% === 4. Loop Through Stimulation Sets & Plot ===
nSets = numel(RespStruct.set);

for sIdx = 1:nSets
    
    % --- IDENTIFY SIGNIFICANT CHANNELS ---
    ptd_idx = 1; 
    
    if amp_idx > numel(RespStruct.set(sIdx).amp)
        warning('Set %d does not have data for Amp Index %d', sIdx, amp_idx);
        continue;
    end
    
    chan_structs = RespStruct.set(sIdx).amp(amp_idx).ptd(ptd_idx).channel;
    
    sig_channels = [];
    for c = 1:length(chan_structs)
        if isfield(chan_structs(c), 'is_responsive') && chan_structs(c).is_responsive
            sig_channels = [sig_channels, c]; 
        end
    end
    
    % --- PLOTTING ---
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.4, 0.3, 0.22, 0.55]);
    hold on;
    
    % --- TITLE LOGIC (MODIFIED) ---
    % 1. Get Stim Channels
    stimCh_List = uniqueComb(sIdx, :);
    stimCh_List = stimCh_List(stimCh_List > 0);
    stimCh_Str  = num2str(stimCh_List);

    % 2. Get Set Name
    if contains(lower(data_folder), 'sim')
        stimName = 'Simultaneous';
    elseif contains(lower(data_folder), 'seq')
        stimName = sprintf('Sequential Set %d', sIdx);
    else
        stimName = sprintf('Single Set %d', sIdx);
    end
    
    % 3. Print Title
    title(sprintf('%s (Ch: %s)\n%.1f µA', stimName, stimCh_Str, amp_to_plot), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    axis equal; axis off;
    label_offset = 1.2;
    
    % Draw Channels
    for s_shank = 1:nShank
        for ch = channel_groups{s_shank}
            idx = find(ch == [channel_groups{:}]);
            if isempty(idx), continue; end
            
            xx = x_pos(idx);
            yy = y_pos(idx);
            
            if ismember(ch, sig_channels)
                markerColor = 'r';      
                edgeColor = 'k';
            else
                markerColor = [0.8 0.8 0.8]; 
                edgeColor = [0.6 0.6 0.6];
            end
            
            plot(xx, yy, 'o', 'MarkerSize', 10, ...
                 'MarkerFaceColor', markerColor, ...
                 'MarkerEdgeColor', edgeColor, 'LineWidth', 1);
             
            text(xx + label_offset, yy, sprintf('%d', ch), ...
                 'FontSize', 8, 'Color', [0.2 0.2 0.2], ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
    end
    
    % Draw Triangular Border
    for s_shank = 1:nShank
        nCh_shank = numel(channel_groups{s_shank});
        tip_y = -2;
        base_y = (nCh_shank + 1) * spacing_y;
        half_width_base = 1.4;
        half_width_tip = 0.2;
        x_center = (s_shank-1) * spacing_x;
        
        x_border = [x_center-half_width_base, x_center+half_width_base, ...
                    x_center+half_width_tip, x_center-half_width_tip];
        y_border = [base_y, base_y, tip_y, tip_y];
        
        patch(x_border, y_border, [0.9 0.9 0.9], ...
              'EdgeColor', [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'LineWidth', 1, ...
              'HandleVisibility','off');
    end
    
    xlim([min(x_pos)-3, max(x_pos)+4]);
    ylim([-3, max(y_pos)+2]);
    hold off;
    
    fprintf('Set %d (Stim: %s): %d significant channels found.\n', sIdx, stimCh_Str, length(sig_channels));
end