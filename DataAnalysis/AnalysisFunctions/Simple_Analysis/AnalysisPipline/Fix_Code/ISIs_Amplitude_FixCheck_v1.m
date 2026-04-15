%% ========================================================================
%  Verification Script: Plot Saved Spike Counts
%  Loads the modified .mat file and generates the ISI Tuning Curves
%  to visually verify that the amplitudes were swapped correctly.
% ========================================================================
clear;

% --- 1. USER SETTINGS ---
data_folder  = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_ISI_SimSeq2';

% Paste the EXACT path to the file you just swapped here:
results_path = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX018/Result_SpikeCount_FixWin_5_10uA_Xia_ISI_SimSeq2.mat';

% --- 2. LOAD DATA ---
fprintf('Loading saved results...\n');
load(results_path, 'ResultFR');

% Extract Metadata
target_Amps = ResultFR.Metadata.TargetAmps;
target_ISIs = ResultFR.Metadata.TargetISIs;
uniqueComb  = ResultFR.Metadata.Stimulation_Sets;
nSets       = size(uniqueComb, 1);

% Extract Matrices
SpikeCount_sim = ResultFR.SpikeCounts.Sim;
SpikeCount_seq = ResultFR.SpikeCounts.Seq;

% --- 3. LOAD HARDWARE MAP (For internal matrix indexing) ---
% We need to know where the hardware placed the Amps and PTDs to index properly
cd(data_folder);
S_exp = load(dir('*_exp_datafile_*.mat').name, 'StimParams', 'simultaneous_stim');
simN = S_exp.simultaneous_stim;

amps_all = cell2mat(S_exp.StimParams(2:end,16)); 
trialAmps = amps_all(1:simN:end);
[Amps, ~, ~] = unique(trialAmps); Amps(Amps==-1) = 0;

if simN > 1
    PTD_all = cell2mat(S_exp.StimParams(3:simN:end,6)); 
else
    PTD_all = zeros(length(trialAmps),1);
end
[PTDs_ms, ~, ~] = unique(PTD_all/1000);

% Find which matrix indices correspond to our target amplitudes
[~, ai_targets] = ismember(target_Amps, Amps);


% --- 4. GENERATE FIGURES (One per Set) ---
fprintf('Generating Verification Plots...\n');
amp_colors = lines(length(target_Amps)); 

for ss = 1:nSets
    figure('Color','w', 'Position',[100 100 800 600]); 
    hold on;
    
    stimCh = uniqueComb(ss,:); 
    stimCh = stimCh(stimCh>0);
    
    % Loop through Amplitudes
    for a_idx = 1:length(target_Amps)
        current_ai  = ai_targets(a_idx);
        current_amp = target_Amps(a_idx);
        
        y_mean = nan(1, length(target_ISIs));
        
        % Gather the data across ISIs
        for p_idx = 1:length(target_ISIs)
            target_isi = target_ISIs(p_idx);
            
            if target_isi == 0
                data_set = squeeze(SpikeCount_sim(:, current_ai, ss));
            else
                p = find(abs(PTDs_ms - target_isi) < 0.001);
                if isempty(p)
                    data_set = NaN;
                else
                    data_set = squeeze(SpikeCount_seq(:, current_ai, ss, p));
                end
            end
            
            data_set = data_set(~isnan(data_set)); % Clean out invalid channels
            
            if ~isempty(data_set)
                y_mean(p_idx) = mean(data_set);
            end
        end
        
        % Plot the line for this Amplitude
        if any(~isnan(y_mean))
            col = amp_colors(a_idx, :);
            lbl = sprintf('%.1f uA', current_amp); 
            
            plot(target_ISIs, y_mean, '-o', 'Color', col, 'LineWidth', 2, ...
                'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', lbl);
        end
    end
    
    % Formatting
    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold', 'FontSize', 12); 
    ylabel('Spikes Count', 'FontWeight','bold', 'FontSize', 12);
    
    % Marked as Verification to avoid confusion with original plots
    title(sprintf('Set %d (Chns: %s)', ss, num2str(stimCh)), 'FontWeight','bold', 'FontSize', 14);
    xticks(sort(target_ISIs));
    box off;
    
    lgd = legend('Location','best','Box','off'); 
    title(lgd, 'Amplitudes');
end

fprintf('Load Complete.\n');