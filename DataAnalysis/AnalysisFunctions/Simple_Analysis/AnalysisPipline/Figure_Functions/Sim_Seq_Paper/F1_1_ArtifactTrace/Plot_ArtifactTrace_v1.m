%% Raw / Denoised / Filtered Trace Viewer (Averaged Artifact Overlay)
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1';

channels_to_plot = 32;                   % channels to plot (Depth_s index)
% Added multiple amplitudes to overlay on the same plot
amps_to_plot     = [3,5,8];           % amplitudes to overlay (µA)
% Defines the two plots you want: 0 (Simultaneous) and 5 (Sequential)
ptd_to_plot      = [0, 5];               % PTDs (ms)
sets_to_plot     = [];                   % stimulation sets, [] means all

% This now dictates how many trials are averaged together
nTrials_to_plot  = 20;                   
% Tightened window to closely match the artifact reference image
plot_window_ms   = [-1 7];              
y_axis_limits    = [-8000 8000];

Electrode_Type   = 2;                    % 0 rigid, 1 flex, 2 4-shank flex

% 'raw'  = amplifier.dat
% 'dn'   = amplifier_dn_sab.dat   (denoised)
% 'mu'   = <base_name>.mu_sab.dat (filtered / MUA)
trace_type = 'raw';    % raw trace

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
S = load(fDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
n_Trials          = S.n_Trials;
E_MAP             = S.E_MAP;

%% ====================== AMPLITUDES ======================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps,~,ampIdx] = unique(trialAmps(:));
Amps(Amps==-1) = 0;
n_AMP = numel(Amps);

%% ====================== PTD ======================
postTrig_all = cell2mat(StimParams(2:end,6));
postTrig = postTrig_all(2:simultaneous_stim:end);
[PTDs,~,ptdIdx] = unique(postTrig);
n_PTD = numel(PTDs);

%% ====================== SEQUENCE-SENSITIVE SETS ======================
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = idx_all(rr); v = v(v>0);
    stimChPerTrial_all{t} = v(:)';
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb,~,combClass] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);
combClass_win = combClass;

%% ====================== APPLY USER FILTERS ======================
if isempty(ptd_to_plot)
    ptd_sel = PTDs;
else
    ptd_sel = ptd_to_plot * 1000;  % user inputs ms, PTDs are in µs
    ptd_sel = intersect(PTDs, ptd_sel);
end

if isempty(sets_to_plot)
    set_sel = 1:nSets;
else
    set_sel = sets_to_plot;
end

if isempty(amps_to_plot)
    amps_sel = Amps;
else
    amps_sel = intersect(Amps, amps_to_plot);
end

%% ====================== TIME SAMPLES ======================
samp_win = round(plot_window_ms/1000 * FS);
Nsamp = samp_win(2)-samp_win(1)+1;
time_ms = (samp_win(1):samp_win(2)) / FS * 1000;

%% ====================== CHANNEL MAP ======================
d = Depth_s(Electrode_Type);

%% ====================== MAIN LOOP ======================
for ich = 1:length(channels_to_plot)

    ch_plot  = channels_to_plot(ich);
    ch_intan = d(ch_plot);

    for i_set = 1:length(set_sel)
        set_id = set_sel(i_set);

        stimVec   = uniqueComb(set_id,:);
        stimVec   = stimVec(stimVec>0);
        set_label = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimVec,'UniformOutput',false),'→');

        % Swapped the loop order: Create Figure per PTD first
        for i_ptd = 1:length(ptd_sel)
            ptd_val = ptd_sel(i_ptd);
            
            if ptd_val == 0
                mode_str = 'Simultaneous';
            else
                mode_str = sprintf('Sequential (%d ms)', ptd_val/1000);
            end

            figName = sprintf('%sTrace_Ch%d_Set%s_%s', data_label, ch_plot, set_label, mode_str);
            % figure('Name',figName,'Color','w','Position',[150 150 800 600]);
            figure('Name',figName,'Color','w','Units','centimeters','Position',[5 5 8.8 8.8]);
            hold on;
            
            % Setup colors for the overlaid amplitudes (Grayscale gradient)
            amp_colors = flipud(gray(length(amps_sel) + 2)); 
            
            legend_handles = gobjects(length(amps_sel), 1);
            legend_labels  = cell(length(amps_sel), 1);

            % Loop through Amplitudes inside the figure
            for i_amp = 1:length(amps_sel)
                amp_val = amps_sel(i_amp);

                tidx_set = find(combClass_win == set_id);
                tidx_ptd = find(ismember(postTrig, ptd_val));
                tidx_amp = find(trialAmps == amp_val);

                tlist = intersect(intersect(tidx_set, tidx_ptd), tidx_amp);

                if isempty(tlist), continue; end
                tlist = tlist(1:min(nTrials_to_plot, numel(tlist)));
                num_valid_trials = numel(tlist);
                
                % Pre-allocate matrix to store trial data for averaging
                all_trial_data = zeros(num_valid_trials, Nsamp);

                for k = 1:num_valid_trials
                    tr = tlist(k);
                    start_idx = trig(tr) + samp_win(1);
                    byte_pos  = start_idx * nChn * 2;

                    fseek(fid, byte_pos, 'bof');
                    data_block = fread(fid, [nChn, Nsamp], 'int16') * 0.195;
                    
                    all_trial_data(k, :) = data_block(ch_intan,:);
                end
                
                % Calculate the mean across the collected trials
                mean_trace = mean(all_trial_data, 1);
                
                % Plot the averaged trace
                col = amp_colors(i_amp + 1, :); 
                legend_handles(i_amp) = plot(time_ms, mean_trace, 'LineWidth', 1.2, 'Color', col);
                legend_labels{i_amp}  = sprintf('%d µA', amp_val);
            end

            % Figure Formatting and dynamic legend
            xline(0,'r-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            % Draw a second red line for the Sequential delay
            if ptd_val > 0
                xline(ptd_val/1000, 'r-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
            % title(sprintf('%s Trace | %s | Rec: Ch %d', data_label, mode_str, ch_plot), 'FontWeight', 'bold', 'FontSize', 14);
            % xlabel('Time (ms)', 'FontWeight', 'bold');
            % ylabel('Voltage (µV)', 'FontWeight', 'bold');

            xlabel('Time (ms)', 'FontSize', 10);
            ylabel('Voltage (µV)', 'FontSize', 10);
            
            % Only show legend for amplitudes that actually existed in the data
            valid_leg = ~arrayfun(@(x) isequal(class(x), 'matlab.graphics.GraphicsPlaceholder'), legend_handles);
            legend(legend_handles(valid_leg), legend_labels(valid_leg), 'Location', 'northeast', 'Box', 'off');
            box off;
            axis square;
            xlim(plot_window_ms);
            ylim(y_axis_limits);
            xticks(-1:1:7); 
            yticks([-8000, -4000, 0, 4000, 8000]);
            set(gca, 'FontSize', 8);
        end
    end
end

fclose(fid);