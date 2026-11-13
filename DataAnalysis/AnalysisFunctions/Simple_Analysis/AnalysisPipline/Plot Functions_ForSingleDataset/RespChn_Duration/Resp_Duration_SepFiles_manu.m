%% ============================================================
%   MANUAL RESPONSE DURATION REVIEW GUI (Draggable ROI)
%   - Loads raw data and previously calculated automatic results.
%   - Presents Raster (top) and Smoothed PSTH (bottom) per channel.
%   - Allows user to adjust Start/End times using a draggable rectangle.
%   - Saves results to a new _ManualReviewed.mat file.
% ============================================================
% NOTE: This is a Script. Do not put 'function' at the top.
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% 1. Raw Data Path (Where the spike files are)
raw_data_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Seq4'; 

% 2. Select the Automatic Result File (Using a Popup)
fprintf('Please select the AUTOMATIC RESULT .mat file from the previous step...\n');
[file_name, file_path] = uigetfile('*.mat', 'Select the Auto-Analysis Result File');

if isequal(file_name, 0)
    error('No file selected. Analysis cancelled.');
end
auto_result_mat = fullfile(file_path, file_name);
fprintf('Selected File: %s\n', auto_result_mat);

% 3. GUI View Settings
view_win_ms = [-10 50];   % X-axis limits for plotting (ms)
psth_bin_ms = 1;          % Bin size for PSTH display
smooth_sigma_ms = 2;      % Gaussian smoothing for PSTH display (easier to see)

% 4. Electrode Type (needed to map channel index to depth)
Electrode_Type = 2;

%% ================= 1. SETUP & LOAD DATA =================
fprintf('Loading Raw Data from: %s\n', raw_data_folder);
[R_raw, sp, trig, S_raw, QC] = load_experiment_data(raw_data_folder);
d_map = Depth_s(Electrode_Type); FS = 30000;

fprintf('Loading Automatic Results...\n');
AutoRes = load(auto_result_mat);
if ~isfield(AutoRes, 'DurationData'), error('Invalid Result MAT file. Does not contain DurationData.'); end

% --- Helper: Re-identify PTD indices from raw data ---
Stim = S_raw.StimParams; simN = S_raw.simultaneous_stim;
if isfield(S_raw, 'n_Trials'), nTr = S_raw.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end
amps_all = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
[~,~,ampIdx] = unique(trialAmps);
if simN > 1, PTD_all_us = cell2mat(Stim(3:simN:end,6)); else, PTD_all_us = zeros(nTr,1); end
PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);
ptd_sim_idx = find(abs(PTDs_ms - 0) < 0.001);
% Assume whatever non-zero PTD exists is the target seq PTD
ptd_seq_idx = find(PTDs_ms > 0.001, 1); 

stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S_raw.E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[~,~,combClass] = unique(comb,'rows','stable'); 

%% ================= 2. FLATTEN DATA FOR NAVIGATION =================
% Convert nested structure into a flat list for easy "Next/Prev" navigation.
fprintf('Flattening data structure for review...\n');
ReviewList = struct('SetIdx',{},'AmpIdx',{},'ModeStr',{},'ChIdx',{},'tStart',{},'tEnd',{});
dd = AutoRes.DurationData;

for ss = 1:length(dd.Set)
    if isempty(dd.Set(ss).Label) || ~isfield(dd.Set(ss), 'Amp'), continue; end
    % Check if .Amp is a struct (Fix for the expansion bug)
    if ~isstruct(dd.Set(ss).Amp), continue; end
    
    for ai = 1:length(dd.Set(ss).Amp)
        % Handle Sim
        if isfield(dd.Set(ss).Amp(ai), 'Sim') && ~isempty(dd.Set(ss).Amp(ai).Sim)
           durs = dd.Set(ss).Amp(ai).Sim;
           resp_chs = [];
           try resp_list = R_raw.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
               for ch=1:length(resp_list), if resp_list(ch).is_responsive, resp_chs=[resp_chs, ch]; end; end
           catch; end
           
           for k = 1:min(length(durs), length(resp_chs))
               idx = length(ReviewList)+1;
               ReviewList(idx).SetIdx = ss; ReviewList(idx).AmpIdx = ai;
               ReviewList(idx).ModeStr = 'Sim'; ReviewList(idx).ChIdx = resp_chs(k);
               ReviewList(idx).tStart = 5; 
               ReviewList(idx).tEnd = 5 + durs(k);
           end
        end
        % Handle Seq
        if isfield(dd.Set(ss).Amp(ai), 'Seq') && ~isempty(dd.Set(ss).Amp(ai).Seq)
           durs = dd.Set(ss).Amp(ai).Seq;
           resp_chs = [];
           try resp_list = R_raw.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
               for ch=1:length(resp_list), if resp_list(ch).is_responsive, resp_chs=[resp_chs, ch]; end; end
           catch; end
           
           for k = 1:min(length(durs), length(resp_chs))
               idx = length(ReviewList)+1;
               ReviewList(idx).SetIdx = ss; ReviewList(idx).AmpIdx = ai;
               ReviewList(idx).ModeStr = 'Seq'; ReviewList(idx).ChIdx = resp_chs(k);
               ReviewList(idx).tStart = 5; 
               ReviewList(idx).tEnd = 5 + durs(k);
           end
        end
    end
end
nTotal = length(ReviewList);
fprintf('Found %d unique channel responses to review.\n', nTotal);
if nTotal == 0, error('No data found to review.'); end

%% ================= 3. Initialize GUI =================
% Create Figure
hFig = figure('Name','Manual Duration Review', 'Position',[100, 100, 1000, 800], 'Color','w');

% Layout
mainLayout = uigridlayout(hFig, [3, 1], 'RowHeight', {'fit', '1fr', 'fit'});

% Top Info Panel
infoPanel = uipanel(mainLayout, 'BorderType','none', 'BackgroundColor','w');
hTitle = uicontrol(infoPanel, 'Style','text', 'String','Initializing...', ...
    'FontSize',14, 'FontWeight','bold', 'Units','normalized', 'Position',[0 0 1 1], 'BackgroundColor','w');

% Plotting Axes
plotLayout = tiledlayout(mainLayout, 2, 1, 'TileSpacing','compact');
axRaster = nexttile(plotLayout); title('Raster Plot');
axPSTH = nexttile(plotLayout); title('Smoothed PSTH'); xlabel('Time (ms)');
linkaxes([axRaster, axPSTH], 'x');

% Bottom Control Panel
cPanel = uigridlayout(mainLayout, [1, 4], 'ColumnWidth', {'1fr', '1fr', '2fr', '1fr'});
btnPrev = uicontrol(cPanel, 'Style','pushbutton', 'String','< Previous', 'FontSize',12);
txtProgress = uicontrol(cPanel, 'Style','text', 'String','1 / N', 'FontSize',12);
btnNext = uicontrol(cPanel, 'Style','pushbutton', 'String','Next >', 'FontSize',12);
btnSave = uicontrol(cPanel, 'Style','pushbutton', 'String','Save Review & Quit', ...
    'FontSize',12, 'FontWeight','bold', 'BackgroundColor',[0.2 0.8 0.2]);

% --- Bundle Data for Callbacks ---
guiData.currIdx = 1;
guiData.ReviewList = ReviewList;
guiData.AutoRes = AutoRes; 
guiData.R_raw = R_raw;
guiData.sp = sp; 
guiData.trig = trig;
guiData.ampIdx = ampIdx;
guiData.combClass = combClass;
guiData.ptdIdx = ptdIdx;
guiData.ptd_sim_idx = ptd_sim_idx;
guiData.ptd_seq_idx = ptd_seq_idx;
guiData.d_map = d_map;
guiData.FS = FS;
guiData.view_win_ms = view_win_ms;
guiData.psth_bin_ms = psth_bin_ms;
guiData.smooth_sigma_ms = smooth_sigma_ms;
guiData.auto_result_mat = auto_result_mat;
guiData.ROI = [];
guiData.isModified = false;

% Handles
guiData.hTitle = hTitle;
guiData.axRaster = axRaster;
guiData.axPSTH = axPSTH;
guiData.txtProgress = txtProgress;
guiData.hFig = hFig;

% Assign Callbacks
set(btnPrev, 'Callback', {@cb_nav, -1});
set(btnNext, 'Callback', {@cb_nav, 1});
set(btnSave, 'Callback', @cb_save);
set(hFig,    'CloseRequestFcn', @cb_close);

% Store Data
set(hFig, 'UserData', guiData);

% Trigger Initial Display
update_display(hFig);


%% ================= LOCAL FUNCTIONS =================
% These functions must be at the END of the file.

function cb_nav(src, ~, delta)
    hFig = ancestor(src, 'figure');
    guiData = get(hFig, 'UserData');
    
    % Save current ROI
    guiData = save_current_roi(guiData);
    
    % Update Index
    newIdx = guiData.currIdx + delta;
    if newIdx < 1 || newIdx > length(guiData.ReviewList), return; end
    
    guiData.currIdx = newIdx;
    set(hFig, 'UserData', guiData);
    update_display(hFig);
end

function cb_save(src, ~)
    hFig = ancestor(src, 'figure');
    guiData = get(hFig, 'UserData');
    guiData = save_current_roi(guiData); % Save last view
    
    ReviewList = guiData.ReviewList;
    UpdatedDurationData = guiData.AutoRes.DurationData;
    R_raw = guiData.R_raw;
    
    fprintf('Reconstructing and saving data...\n');
    nTotal = length(ReviewList);
    
    for i = 1:nTotal
        item = ReviewList(i);
        duration = item.tEnd - item.tStart;
        
        if strcmp(item.ModeStr, 'Sim')
            % Find channel index in the list
            try 
                resp_list = R_raw.set(item.SetIdx).amp(item.AmpIdx).ptd(guiData.ptd_sim_idx).channel;
                resp_chs = [];
                for ch=1:length(resp_list), if resp_list(ch).is_responsive, resp_chs=[resp_chs, ch]; end; end
                idx_in_list = find(resp_chs == item.ChIdx);
                if ~isempty(idx_in_list)
                    UpdatedDurationData.Set(item.SetIdx).Amp(item.AmpIdx).Sim(idx_in_list) = duration;
                end
            catch; end
        else
            try
                resp_list = R_raw.set(item.SetIdx).amp(item.AmpIdx).ptd(guiData.ptd_seq_idx).channel;
                resp_chs = [];
                for ch=1:length(resp_list), if resp_list(ch).is_responsive, resp_chs=[resp_chs, ch]; end; end
                idx_in_list = find(resp_chs == item.ChIdx);
                if ~isempty(idx_in_list)
                    UpdatedDurationData.Set(item.SetIdx).Amp(item.AmpIdx).Seq(idx_in_list) = duration;
                end
            catch; end
        end
    end
    
    % Save
    [fPath, fName, fExt] = fileparts(guiData.auto_result_mat);
    saveName = fullfile(fPath, [fName '_ManualReviewed' fExt]);
    
    % Handle saving structure safely
    AutoRes = guiData.AutoRes;
    AutoRes.DurationData = UpdatedDurationData;
    save(saveName, '-struct', 'AutoRes');
    
    fprintf('>>> Manual Review Saved to: %s\n', saveName);
    delete(hFig);
end

function cb_close(src, ~)
    guiData = get(src, 'UserData');
    if guiData.isModified
       choice = questdlg('You have unsaved changes. Quit anyway?', 'Unsaved Changes', 'Quit (Discard)', 'Cancel', 'Cancel');
       if strcmp(choice, 'Quit (Discard)')
           delete(src);
       end
    else
        delete(src);
    end
end

function guiData = save_current_roi(guiData)
    if ~isempty(guiData.ROI) && isvalid(guiData.ROI)
        pos = guiData.ROI.Position;
        tStart = pos(1);
        tEnd = pos(1) + pos(3);
        guiData.ReviewList(guiData.currIdx).tStart = tStart;
        guiData.ReviewList(guiData.currIdx).tEnd = tEnd;
        guiData.isModified = true;
    end
end

function update_display(hFig)
    guiData = get(hFig, 'UserData');
    item = guiData.ReviewList(guiData.currIdx);
    
    % Update Title
    set(guiData.txtProgress, 'String', sprintf('%d / %d', guiData.currIdx, length(guiData.ReviewList)));
    ampVal = guiData.AutoRes.DurationData.Set(item.SetIdx).Amp(item.AmpIdx).Val;
    set(guiData.hTitle, 'String', sprintf('Set %d | %.1f uA | %s | Ch %d (Depth %d)', ...
        item.SetIdx, ampVal, item.ModeStr, item.ChIdx, guiData.d_map(item.ChIdx)));
    
    % Get Trials
    if strcmp(item.ModeStr, 'Sim')
        pidx = guiData.ptd_sim_idx;
    else
        pidx = guiData.ptd_seq_idx;
    end
    tr_ids = find(guiData.ampIdx == item.AmpIdx & guiData.combClass == item.SetIdx & guiData.ptdIdx == pidx);
    
    % Calc PSTH
    edges = guiData.view_win_ms(1) : guiData.psth_bin_ms : guiData.view_win_ms(2);
    centers = edges(1:end-1) + guiData.psth_bin_ms/2;
    raster_x = []; raster_y = []; all_spikes = [];
    
    sp_data = guiData.sp{guiData.d_map(item.ChIdx)};
    
    for k = 1:length(tr_ids)
        tr = tr_ids(k); t0 = guiData.trig(tr)/guiData.FS*1000;
        win_pad = [guiData.view_win_ms(1)-5, guiData.view_win_ms(2)+5];
        tt = sp_data(:,1) - t0;
        mask = tt >= win_pad(1) & tt <= win_pad(2);
        sp_trial = tt(mask);
        all_spikes = [all_spikes; sp_trial]; %#ok<AGROW>
        raster_x = [raster_x; sp_trial]; %#ok<AGROW>
        raster_y = [raster_y; repmat(k, length(sp_trial), 1)]; %#ok<AGROW>
    end
    
    counts = histcounts(all_spikes, edges);
    if guiData.smooth_sigma_ms > 0
        k_w = 3*guiData.smooth_sigma_ms; 
        gk = exp(-(-k_w:guiData.psth_bin_ms:k_w).^2/(2*guiData.smooth_sigma_ms^2)); 
        gk = gk/sum(gk);
        smooth_counts = conv(counts, gk, 'same');
    else
        smooth_counts = counts;
    end
    
    % Plot Raster
    ax = guiData.axRaster;
    cla(ax); axes(ax); hold on;
    plot(raster_x, raster_y, 'k.', 'MarkerSize', 4);
    xlim(guiData.view_win_ms); ylim([0 length(tr_ids)+1]);
    xline(0, 'r--');
    
    % Plot PSTH
    ax = guiData.axPSTH;
    cla(ax); axes(ax); hold on;
    plot(centers, smooth_counts, 'k-', 'LineWidth', 1.5);
    xlim(guiData.view_win_ms); grid on;
    xline(0, 'r--');
    
    % Draw ROI
    roi_x = item.tStart;
    roi_w = item.tEnd - item.tStart;
    yl = ylim(ax);
    
    roi = drawrectangle(ax, 'Position', [roi_x, yl(1), roi_w, yl(2)-yl(1)], ...
        'Color', [0 0.4 0.8], 'FaceAlpha', 0.2, 'LineWidth', 1.5, ...
        'Label', sprintf('%.1f ms', roi_w), 'LabelVisible', 'on');
    
    addlistener(roi, 'MovingROI', @roi_moving);
    guiData.ROI = roi;
    set(guiData.hFig, 'UserData', guiData);
end

function roi_moving(src, evt)
    pos = evt.CurrentPosition;
    src.Label = sprintf('%.1f ms', pos(3));
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*.BadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end