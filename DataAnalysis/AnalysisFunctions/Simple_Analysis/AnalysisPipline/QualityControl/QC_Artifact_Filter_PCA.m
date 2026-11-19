%% ============================================================
% SpikeFiltering_PCA_Cleanup.m
% Refined: PCA-Based Outlier Rejection (Shape Analysis)
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Sim';
FS = 30000;                 

% 1. FILTERING PARAMETERS
PCA_Outlier_Sigma    = 5;     % Keep spikes within N standard deviations of baseline cluster
Min_Baseline_Spikes  = 10;    % Minimum spikes needed to learn the shape
Max_Biological_uV    = 500;   % Safety net (just in case of saturation)

% 2. WINDOWS
baseline_window_ms   = [-50 -5]; % The "Safe Zone" (Learns real spike shape)
late_evoked_window   = [10 50];   % Fallback safe zone
plot_diagnostics     = 1;         % 1 = Show plots for every channel

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder), error('Folder does not exist.'); end
cd(data_folder);
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4, base_name = last_folder(1 : underscores(end-1)-1);
else, base_name = last_folder; end

%% ================= LOAD SPIKES =================
fname_sp = [base_name '.sp_xia_FirstPulse.mat'];
assert(isfile(fname_sp), 'Cannot find %s.', fname_sp);
S_in = load(fname_sp);
if isfield(S_in,'sp_seq'), sp_in = S_in.sp_seq;
elseif isfield(S_in,'sp'), sp_in = S_in.sp;
else, error('No sp or sp_clipped found in %s', fname_sp); end

nCh = numel(sp_in);
if isempty(dir('*.trig.dat')), cur=pwd; cleanTrig_sabquick; cd(cur); end
trig = loadTrig(0);
nTrials = numel(trig);

%% ================= PCA OUTLIER REJECTION & PLOTTING =================
sp_pca = sp_in; % Output structure
fprintf('\nRunning PCA-Based Artifact Rejection (Sigma = %.1f)...\n', PCA_Outlier_Sigma);

for ch = 1:nCh
    if isempty(sp_in{ch}), continue; end
    spt = sp_in{ch}(:,1);
    wfs = sp_in{ch}(:,2:end);
    
    nSpikes = size(wfs, 1);
    if nSpikes < 20
        fprintf('Ch %2d: Too few spikes (%d) for PCA. Skipping.\n', ch, nSpikes);
        continue; 
    end
    
    % --- 1. Identify Training Set (Baseline Spikes) ---
    mask_base = false(nSpikes, 1);
    for tr = 1:nTrials
        t0 = trig(tr)/FS*1000;
        mask_base = mask_base | (spt >= t0 + baseline_window_ms(1) & spt <= t0 + baseline_window_ms(2));
    end
    
    % Fallback to late evoked if baseline is empty
    if sum(mask_base) < Min_Baseline_Spikes
        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            mask_base = mask_base | (spt >= t0 + late_evoked_window(1) & spt <= t0 + late_evoked_window(2));
        end
    end
    
    if sum(mask_base) < Min_Baseline_Spikes
        fprintf('Ch %2d: Not enough baseline spikes to learn shape. Skipping.\n', ch);
        continue;
    end
    
    % --- 2. Run PCA ---
    % Project all waveforms into 2D shape space
    [coeff, score, ~] = pca(wfs, 'NumComponents', 2);
    
    % --- 3. Define "Real Spike" Cluster ---
    % Get stats of the baseline spikes in PC space
    score_base = score(mask_base, :);
    mu         = median(score_base); % Use median for robustness
    sigma      = cov(score_base);
    
    % Calculate Mahalanobis Distance for ALL spikes from the baseline center
    % (Distance normalized by the cluster's spread)
    d_sq = mahalanobis_dist(score, mu, sigma);
    dist_sigma = sqrt(d_sq);
    
    % --- 4. Filter ---
    % Reject if distance > Sigma Threshold OR Amplitude > Safety Limit
    peak_vals  = max(abs(wfs), [], 2);
    is_outlier = (dist_sigma > PCA_Outlier_Sigma) | (peak_vals > Max_Biological_uV);
    
    % Save clean spikes
    sp_pca{ch} = sp_in{ch}(~is_outlier, :);
    removed_count = sum(is_outlier);
    
    fprintf('Ch %2d: Removed %d artifacts (%.1f%%).\n', ch, removed_count, (removed_count/nSpikes)*100);
    
    % --- 5. DIAGNOSTIC PLOT ---
    if plot_diagnostics && removed_count > 0
        figure('Color','w', 'Name', sprintf('Ch %d PCA Cleaning', ch), 'Position', [100 100 1200 400]);
        tiledlayout(1,4, 'Padding', 'compact');
        
        % Plot A: PCA Cluster
        nexttile; hold on;
        scatter(score(~is_outlier,1), score(~is_outlier,2), 10, 'b', 'filled', 'MarkerFaceAlpha',0.3);
        scatter(score(is_outlier,1), score(is_outlier,2), 20, 'r', 'filled');
        xlabel('PC 1'); ylabel('PC 2'); title('PCA Shape Space');
        legend('Real', 'Artifact', 'Location','best'); grid on; axis square;
        
        % Plot B: Histogram of Timing (Check if artifacts are at 6-8ms)
        nexttile; hold on;
        histogram(spt(~is_outlier, 1), 50, 'FaceColor','b', 'EdgeColor','none');
        histogram(spt(is_outlier, 1), 50, 'FaceColor','r', 'EdgeColor','none');
        xlabel('Time (ms)'); title('Spike Timing');
        
        % Plot C: Artifact Waveforms
        nexttile; 
        plot(wfs(is_outlier,:)', 'r', 'Color', [1 0 0 0.1]);
        title(sprintf('Artifacts (%d)', removed_count)); ylim([-800 800]); xlim([0 48]);
        
        % Plot D: Clean Waveforms
        nexttile;
        plot(wfs(~is_outlier,:)', 'b', 'Color', [0 0 1 0.05]);
        title(sprintf('Clean (%d)', sum(~is_outlier))); ylim([-300 300]); xlim([0 48]);
        
        drawnow;
    end
end

%% ================= SAVE OUTPUT =================
QC_params = struct();
QC_params.Method            = 'PCA_Mahalanobis';
QC_params.PCA_Sigma         = PCA_Outlier_Sigma;
QC_params.Baseline_Window   = baseline_window_ms;
QC_params.Max_Biological_uV = Max_Biological_uV;

out_name = [base_name '.sp_xia_SSD.mat']; % Keeping same output name for compatibility
if isfile(out_name), delete(out_name); end 
save(out_name, 'sp_in','sp_pca','QC_params','-v7.3');
fprintf('\nSaved cleaned spikes to %s\n', out_name);

%% ================= HELPER FUNCTION =================
function d2 = mahalanobis_dist(X, mu, sigma)
    % Robust distance calculation
    d = X - mu;
    % If covariance is singular (e.g. flat line), add tiny jitter
    if rcond(sigma) < 1e-10
        sigma = sigma + eye(size(sigma))*1e-6; 
    end
    d2 = sum((d / sigma) .* d, 2);
end