%% ============================================================
%   IMPLANTATION ORIENTATION CONTROL ANALYSIS
%
%   Question:
%   Does the sequential stimulation advantage depend on whether
%   the probe was implanted approximately perpendicular or
%   tangential to the cortical surface?
%
%   Metric:
%       Delta = Sequential - Simultaneous
%
%   Analysis:
%       1. Calculate set-level Sim and Seq mean response at each amplitude
%       2. Calculate set-level Delta = Seq - Sim
%       3. Compare Delta between orientations at each amplitude
%          using Wilcoxon rank-sum test
%       4. Collapse across amplitudes within each stimulation set
%       5. Collapse stimulation sets within each animal
%       6. Compare overall animal-level Delta between orientations
%
%   Assumed ResultNorm structure:
%       Norm_Sim : Channel x Amplitude x Set
%       Norm_Seq : Channel x Amplitude x Set x Order
%                  or Channel x Amplitude x Set
%
% ============================================================

clear;
clc;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================

% ------------------------------------------------------------
% PUT YOUR PERPENDICULAR / VERTICAL IMPLANTATION FILES HERE
% ------------------------------------------------------------
perpendicular_files = {
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
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX005/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
};

% ------------------------------------------------------------
% PUT YOUR TANGENTIAL / PARALLEL IMPLANTATION FILES HERE
% ------------------------------------------------------------
tangential_files = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX018/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX018/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX018/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX018/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_4.mat';
    
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
};

% Amplitudes to analyse
target_amps = [1 2 3 4 5 6 8 10];

% Minimum N in EACH orientation group for amplitude comparison
min_n_per_group = 1;

% Overall summary:
% mean Delta across amplitudes within each stimulation set
summary_method = 'mean';

%% ================= 2. PROCESS BOTH ORIENTATIONS =================

fprintf('\n============================================================\n');
fprintf('PROCESSING IMPLANTATION ORIENTATION DATA\n');
fprintf('============================================================\n');

Perp = process_orientation(perpendicular_files, 1, target_amps);
Tang = process_orientation(tangential_files, 2, target_amps);

% Columns in .SetData:
% 1 = Animal ID
% 2 = File ID
% 3 = Set ID
% 4 = Amplitude
% 5 = Sim response
% 6 = Seq response
% 7 = Delta (Seq - Sim)
% 8 = Orientation
%     1 = perpendicular
%     2 = tangential

AllSetData = [Perp.SetData; Tang.SetData];

fprintf('\nTotal valid set-amplitude observations:\n');
fprintf('  Perpendicular : %d\n', size(Perp.SetData,1));
fprintf('  Tangential    : %d\n', size(Tang.SetData,1));

%% ================= 3. ANIMAL-LEVEL ANALYSIS AT EACH AMPLITUDE =================

fprintf('\n============================================================\n');
fprintf('ANIMAL-LEVEL ORIENTATION COMPARISON BY AMPLITUDE\n');
fprintf('Delta = Sequential - Simultaneous\n');
fprintf('Positive Delta = Sequential response is larger\n');
fprintf('============================================================\n');

% Bonferroni correction across amplitudes
bonferroni_n = length(target_amps);

% Results structure
AmpResults = struct();

fprintf(['%-6s | %-27s | %-27s | %-10s | %-10s\n'], ...
    'Amp', 'Perpendicular', 'Tangential', 'p raw', 'p corr');

fprintf(['-----------------------------------------------------------------------------------------------------\n']);

for a = 1:length(target_amps)

    amp = target_amps(a);

    %% ---------------------------------------------------------
    % 1. Get all set-level data at this amplitude
    %% ---------------------------------------------------------

    perp_amp_data = Perp.SetData(abs(Perp.SetData(:,4) - amp) < 0.001, :);
    tang_amp_data = Tang.SetData(abs(Tang.SetData(:,4) - amp) < 0.001, :);

    %% ---------------------------------------------------------
    % 2. Average stimulation sets WITHIN each animal
    %% ---------------------------------------------------------

    perp_animals = animal_average_at_amp(perp_amp_data);
    tang_animals = animal_average_at_amp(tang_amp_data);

    % Columns:
    % 1 = Animal ID
    % 2 = Mean Delta for that animal
    % 3 = Number of stimulation sets contributing

    delta_perp = perp_animals(:,2);
    delta_tang = tang_animals(:,2);

    %% ---------------------------------------------------------
    % 3. Descriptive statistics across animals
    %% ---------------------------------------------------------

    mean_perp = mean(delta_perp, 'omitnan');
    sem_perp  = sem_calc(delta_perp);
    n_perp    = length(delta_perp);

    mean_tang = mean(delta_tang, 'omitnan');
    sem_tang  = sem_calc(delta_tang);
    n_tang    = length(delta_tang);

    %% ---------------------------------------------------------
    % 4. Compare implantation orientations
    %% ---------------------------------------------------------

    if n_perp >= min_n_per_group && n_tang >= min_n_per_group

        % Independent groups:
        % perpendicular animals vs tangential animals
        p_raw = ranksum(delta_perp, delta_tang);

        % Bonferroni correction
        p_corr = min(p_raw * bonferroni_n, 1);

    else
        p_raw  = NaN;
        p_corr = NaN;
    end

    %% ---------------------------------------------------------
    % 5. Optional: test whether Delta > 0 within each orientation
    %% ---------------------------------------------------------

    if n_perp >= 2
        p_perp_vs0 = signrank(delta_perp, 0);
    else
        p_perp_vs0 = NaN;
    end

    if n_tang >= 2
        p_tang_vs0 = signrank(delta_tang, 0);
    else
        p_tang_vs0 = NaN;
    end

    %% ---------------------------------------------------------
    % 6. Save results
    %% ---------------------------------------------------------

    AmpResults(a).Amp = amp;

    AmpResults(a).Perp_Mean = mean_perp;
    AmpResults(a).Perp_SEM  = sem_perp;
    AmpResults(a).Perp_N    = n_perp;
    AmpResults(a).Perp_Data = perp_animals;

    AmpResults(a).Tang_Mean = mean_tang;
    AmpResults(a).Tang_SEM  = sem_tang;
    AmpResults(a).Tang_N    = n_tang;
    AmpResults(a).Tang_Data = tang_animals;

    AmpResults(a).P_Raw  = p_raw;
    AmpResults(a).P_Corr = p_corr;

    AmpResults(a).P_Perp_vs0 = p_perp_vs0;
    AmpResults(a).P_Tang_vs0 = p_tang_vs0;

    %% ---------------------------------------------------------
    % 7. Print summary
    %% ---------------------------------------------------------

    fprintf('%4.1f   | %7.4f +/- %-7.4f (N=%-2d) | %7.4f +/- %-7.4f (N=%-2d) | %-10s | %-10s\n', ...
        amp, ...
        mean_perp, sem_perp, n_perp, ...
        mean_tang, sem_tang, n_tang, ...
        format_p(p_raw), format_p(p_corr));

end


%% ================= 4. DETAILED ANIMAL VALUES =================

fprintf('\n\n============================================================\n');
fprintf('ANIMAL-LEVEL DETAILS BY AMPLITUDE\n');
fprintf('============================================================\n');

for a = 1:length(target_amps)

    amp = target_amps(a);

    fprintf('\n---------------- %.1f uA ----------------\n', amp);

    fprintf('\nPerpendicular:\n');

    P = AmpResults(a).Perp_Data;

    for i = 1:size(P,1)

        fprintf('DX%-3d | Delta = %.4f | N sets = %d\n', ...
            P(i,1), P(i,2), P(i,3));

    end


    fprintf('\nTangential:\n');

    T = AmpResults(a).Tang_Data;

    for i = 1:size(T,1)

        fprintf('DX%-3d | Delta = %.4f | N sets = %d\n', ...
            T(i,1), T(i,2), T(i,3));

    end

end


%% ================= 5. WITHIN-ORIENTATION CHECK =================

fprintf('\n\n============================================================\n');
fprintf('SEQUENTIAL ADVANTAGE WITHIN EACH ORIENTATION\n');
fprintf('Animal-level Delta tested against zero\n');
fprintf('============================================================\n');

fprintf('%-6s | %-15s | %-15s\n', ...
    'Amp', 'Perpendicular', 'Tangential');

fprintf('----------------------------------------------------------\n');

for a = 1:length(target_amps)

    fprintf('%4.1f   | %-15s | %-15s\n', ...
        target_amps(a), ...
        format_p(AmpResults(a).P_Perp_vs0), ...
        format_p(AmpResults(a).P_Tang_vs0));

end


%% ================= HELPER FUNCTION =================

function AnimalData = animal_average_at_amp(SetData)

    % Input SetData columns:
    % 1 = Animal ID
    % 2 = File ID
    % 3 = Set ID
    % 4 = Amplitude
    % 5 = Sim
    % 6 = Seq
    % 7 = Delta
    % 8 = Orientation
    %
    % Output:
    % 1 = Animal ID
    % 2 = Mean Delta across stimulation sets
    % 3 = Number of stimulation sets

    AnimalData = [];

    if isempty(SetData)
        return;
    end

    animal_ids = unique(SetData(:,1));

    for i = 1:length(animal_ids)

        animal = animal_ids(i);

        idx = SetData(:,1) == animal;

        delta_values = SetData(idx,7);
        delta_values = delta_values(~isnan(delta_values));

        if isempty(delta_values)
            continue;
        end

        animal_delta = mean(delta_values, 'omitnan');

        AnimalData = [AnimalData; ...
            animal, animal_delta, length(delta_values)];

    end
end


function s = sem_calc(x)

    x = x(~isnan(x));

    if isempty(x)
        s = NaN;

    elseif length(x) == 1
        s = NaN;

    else
        s = std(x) / sqrt(length(x));
    end
end


function p_str = format_p(p)

    if isnan(p)

        p_str = 'NaN';

    elseif p < 0.0001

        p_str = sprintf('%.12e', p);

    else

        p_str = sprintf('%.14f', p);

    end
end

function Out = process_orientation(file_paths, orientation_code, target_amps)

    Out.SetData = [];

    % Columns:
    % 1 Animal ID
    % 2 File ID
    % 3 Set ID
    % 4 Amp
    % 5 Sim
    % 6 Seq
    % 7 Delta
    % 8 Orientation

    for f = 1:length(file_paths)

        if ~exist(file_paths{f}, 'file')
            fprintf('Missing file: %s\n', file_paths{f});
            continue;
        end

        D = load(file_paths{f});

        if ~isfield(D,'ResultNorm')
            fprintf('Skipped: no ResultNorm in %s\n', file_paths{f});
            continue;
        end

        R = D.ResultNorm;

        if ~isfield(R,'Amps') || ...
           ~isfield(R,'Norm_Sim') || ...
           ~isfield(R,'Norm_Seq')
            fprintf('Skipped: missing required fields in %s\n', file_paths{f});
            continue;
        end

        % Extract animal number from DX###
        tok = regexp(file_paths{f}, 'DX(\d+)', 'tokens', 'once');

        if isempty(tok)
            animal_id = f;
        else
            animal_id = str2double(tok{1});
        end

        amps = R.Amps(:);

        sim = R.Norm_Sim;
        seq = R.Norm_Seq;

        % ------------------------------------------------------
        % Expected:
        % Sim = Channel x Amp x Set
        % Seq = Channel x Amp x Set x Order
        %
        % Average recording channels first.
        % For Seq, also average stimulation order.
        % ------------------------------------------------------

        sim_size = size(sim);
        seq_size = size(seq);

        if ndims(sim) < 3
            sim = reshape(sim, size(sim,1), size(sim,2), 1);
        end

        if ndims(seq) < 3
            seq = reshape(seq, size(seq,1), size(seq,2), 1);
        end

        n_sets_sim = size(sim,3);
        n_sets_seq = size(seq,3);
        n_sets = min(n_sets_sim, n_sets_seq);

        for ss = 1:n_sets

            for a = 1:length(target_amps)

                amp = target_amps(a);

                amp_idx = find(abs(amps - amp) < 0.001, 1);

                if isempty(amp_idx)
                    continue;
                end

                % --------------------------
                % Simultaneous
                % --------------------------
                sim_vals = squeeze(sim(:,amp_idx,ss));
                sim_mean = mean(sim_vals(:), 'omitnan');

                % --------------------------
                % Sequential
                % --------------------------
                if ndims(seq) >= 4
                    seq_vals = squeeze(seq(:,amp_idx,ss,:));
                else
                    seq_vals = squeeze(seq(:,amp_idx,ss));
                end

                seq_mean = mean(seq_vals(:), 'omitnan');

                if isnan(sim_mean) || isnan(seq_mean)
                    continue;
                end

                delta = seq_mean - sim_mean;

                Out.SetData = [Out.SetData; ...
                    animal_id, f, ss, amp, ...
                    sim_mean, seq_mean, delta, orientation_code];
            end
        end
    end
end
