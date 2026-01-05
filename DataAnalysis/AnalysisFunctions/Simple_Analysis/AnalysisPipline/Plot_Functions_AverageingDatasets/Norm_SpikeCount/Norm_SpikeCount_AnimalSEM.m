%% ============================================================
%   GRAND AVERAGE: NORMALIZED RECRUITMENT CURVE (N = Animals)
%   - Logic: 
%       1. Group Datasets by Animal ID (e.g., DX015).
%       2. Average all datasets within an animal -> Animal Mean.
%       3. Average across Animals -> Grand Mean +/- SEM.
%   - Result: Valid N (e.g., N=9) instead of pseudoreplication (N=47).
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX005/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
};

save_figure = true;
save_dir    = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/';
fig_name    = 'GrandAverage_Recruitment_Curve_ByAnimal.fig';

%% ================= 2. AGGREGATE BY ANIMAL =================
fprintf('Processing %d datasets...\n', length(file_paths));

% We use a struct to hold data for each animal
% Animals.(ID).Sim = [Amp, Val; Amp, Val; ...]
Animals = struct();

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    
    % --- Extract Animal ID ---
    % Looks for pattern 'DX' followed by digits
    tokens = regexp(file_paths{i}, '(DX\d+)', 'tokens');
    if isempty(tokens)
        warning('Could not parse Animal ID from: %s', file_paths{i});
        animal_id = 'Unknown';
    else
        animal_id = tokens{1}{1};
    end
    
    % Ensure struct field exists
    if ~isfield(Animals, animal_id)
        Animals.(animal_id).Sim = [];
        Animals.(animal_id).Seq = [];
    end
    
    % --- Load Data ---
    D = load(file_paths{i});
    if ~isfield(D, 'ResultNorm'), continue; end
    R = D.ResultNorm;
    
    Amps = R.Amps;
    
    % Collapse Channels & Sets to get Mean for this specific Dataset
    % Sim
    sim_d_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    % Seq (Handle PTD dim if present)
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_d_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_d_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end
    
    % --- Store in Animal's List ---
    for a = 1:length(Amps)
        % Sim
        if ~isnan(sim_d_mean(a))
            Animals.(animal_id).Sim = [Animals.(animal_id).Sim; Amps(a), sim_d_mean(a)];
        end
        % Seq
        if ~isnan(seq_d_mean(a))
            Animals.(animal_id).Seq = [Animals.(animal_id).Seq; Amps(a), seq_d_mean(a)];
        end
    end
end

%% ================= 3. CALCULATE ANIMAL MEANS =================
% Now we collapse the lists inside each animal to get 1 curve per animal.
% Then we add that single curve to the Grand Pool.

Pool_Sim = []; 
Pool_Seq = [];

animal_names = fieldnames(Animals);
fprintf('Found %d Animals: %s\n', length(animal_names), strjoin(animal_names, ', '));

for i = 1:length(animal_names)
    id = animal_names{i};
    
    % --- Sim ---
    data = Animals.(id).Sim;
    if ~isempty(data)
        u_amps = unique(data(:,1));
        for a = 1:length(u_amps)
            amp = u_amps(a);
            % Average all datasets for this animal at this amplitude
            mu = mean(data(data(:,1)==amp, 2)); 
            Pool_Sim = [Pool_Sim; amp, mu]; %#ok<*AGROW>
        end
    end
    
    % --- Seq ---
    data = Animals.(id).Seq;
    if ~isempty(data)
        u_amps = unique(data(:,1));
        for a = 1:length(u_amps)
            amp = u_amps(a);
            mu = mean(data(data(:,1)==amp, 2));
            Pool_Seq = [Pool_Seq; amp, mu];
        end
    end
end

%% ================= 4. GRAND STATISTICS (N = Animals) =================
Unique_Amps = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Amps = sort(Unique_Amps);

Grand_Sim_Mean = nan(size(Unique_Amps)); Grand_Sim_SEM = nan(size(Unique_Amps));
Grand_Seq_Mean = nan(size(Unique_Amps)); Grand_Seq_SEM = nan(size(Unique_Amps));
N_Sim = zeros(size(Unique_Amps)); N_Seq = zeros(size(Unique_Amps));

fprintf('\n=== GRAND STATISTICS (N = Animals) ===\n');
fprintf('Amp(uA)\tN_Sim\tMean_Sim\tSEM_Sim\t\tN_Seq\tMean_Seq\tSEM_Seq\n');

for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    
    % --- Sim ---
    vals = Pool_Sim(Pool_Sim(:,1) == amp, 2);
    n = length(vals);
    if n > 0
        Grand_Sim_Mean(k) = mean(vals);
        % Grand_Sim_SEM(k)  = std(vals) / sqrt(n); % N is now number of animals
        Grand_Sim_SEM(k)  = std(vals); % N is now number of animals
        N_Sim(k) = n;
    end
    
    % --- Seq ---
    vals = Pool_Seq(Pool_Seq(:,1) == amp, 2);
    n = length(vals);
    if n > 0
        Grand_Seq_Mean(k) = mean(vals);
        % Grand_Seq_SEM(k)  = std(vals) / sqrt(n);
        Grand_Seq_SEM(k)  = std(vals);
        N_Seq(k) = n;
    end
    
    fprintf('%.1f\t%d\t%.4f\t%.4f\t\t%d\t%.4f\t%.4f\n', ...
        amp, N_Sim(k), Grand_Sim_Mean(k), Grand_Sim_SEM(k), N_Seq(k), Grand_Seq_Mean(k), Grand_Seq_SEM(k));
end

%% ================= 5. PLOT =================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% Helper Plot
plot_shaded(Unique_Amps, Grand_Sim_Mean, Grand_Sim_SEM, [0 0.4 0.8]);
plot_shaded(Unique_Amps, Grand_Seq_Mean, Grand_Seq_SEM, [0.8 0 0.2]);

% Lines
valid = ~isnan(Grand_Sim_Mean);
plot(Unique_Amps(valid), Grand_Sim_Mean(valid), '-o', 'Color', [0 0.3 0.8], 'LineWidth', 2, 'MarkerFaceColor', [0 0.3 0.8], 'DisplayName', 'Simultaneous');

valid = ~isnan(Grand_Seq_Mean);
plot(Unique_Amps(valid), Grand_Seq_Mean(valid), '-s', 'Color', [0.7 0 0], 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

yline(1.0, '--k', 'Ref (5uA)', 'HandleVisibility','off');
xlabel('Amplitude (\muA)', 'FontSize', 12, 'FontWeight','bold');
ylabel('Normalized Spike Count', 'FontSize', 12, 'FontWeight','bold');
title(sprintf('Grand Average (N=%d Animals)', length(animal_names)), 'FontSize', 14);
legend('Location','best','Box','off');
box off; set(gca, 'FontSize', 12);
xlim([min(Unique_Amps)-0.5, max(Unique_Amps)+0.5]);

% Save
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(gcf, fullfile(save_dir, fig_name));
    fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
end

function plot_shaded(x, y, err, col)
    valid = ~isnan(y) & ~isnan(err);
    if sum(valid) < 2, return; end
    x = x(valid); y = y(valid); err = err(valid);
    fill([x; flipud(x)], [y+err; flipud(y-err)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end