%% ============================================================
%   Raster + PSTH for MULTIPLE SEQUENTIAL PTDs
%   One figure per stimulation set (content-matched across PTDs)
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
seq_folders = {
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_8ms_251125_113541'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_10ms_251125_114346'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_15ms_251125_115256'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_20ms_251125_120014'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_25ms_251125_120717'
};

Electrode_Type = 1;
target_channel = 18;
plot_amps      = [10];

ras_win        = [-10 50];
bin_ms         = 1;
smooth_ms      = 2;
FS             = 30000;

%% Colors
cmapPTD = lines(numel(seq_folders));

%% ===================== PSTH KERNEL ===================== %%
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2); g=g/sum(g);

% LOAD ALL PTD FOLDERS
Seq = struct;
allStimSetStrings = {};  % to collect all unique stimulation sets

for f = 1:numel(seq_folders)
    folder = seq_folders{f};
    cd(folder);
    % Spike data
    if ~isempty(dir('*sp_xia_FirstPulse.mat'))
        load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq');
        Seq(f).sp = sp_seq;
    else
        load(dir('*sp_xia.mat').name,'sp_clipped');
        Seq(f).sp = sp_clipped;
    end
    % Trigger
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
    Seq(f).trig = loadTrig(0);
    % StimParams
    S = load(dir('*_exp_datafile_*.mat').name, ...
             'StimParams','simultaneous_stim','E_MAP','n_Trials');
    Stim = S.StimParams;
    simN = S.simultaneous_stim;
    Seq(f).E_MAP    = S.E_MAP;
    Seq(f).nTrials  = S.n_Trials;
    Seq(f).simN     = simN;
    % Amplitude
    amps_all = cell2mat(Stim(2:end,16));
    Seq(f).trialAmps = amps_all(1:simN:end);
    % PTD (us → ms)
    PTD_all = cell2mat(Stim(2:end,6));
    Seq(f).PTD = PTD_all(2:simN:end)/1000;
    % Ordered stimulation sets
    stimNames = Stim(2:end,1);
    [~, idx] = ismember(stimNames, Seq(f).E_MAP(2:end));
    stimChPerTrial = cell(Seq(f).nTrials,1);
    for t = 1:Seq(f).nTrials
        rr = (t-1)*simN + (1:simN);
        v = idx(rr);
        v = v(v>0);
        stimChPerTrial{t} = v(:).';
    end
    Seq(f).stimChPerTrial = stimChPerTrial;
    comb = zeros(Seq(f).nTrials,simN);
    for t = 1:Seq(f).nTrials
        v = stimChPerTrial{t};
        comb(t,1:numel(v)) = v;
    end
    [Seq(f).uniqueComb,~,Seq(f).combClass] = unique(comb,'rows','stable');
    Seq(f).nSets = size(Seq(f).uniqueComb,1);
    % ---- Collect unique "stimset strings" for matching later ----
    for u = 1:Seq(f).nSets
        setvec = Seq(f).uniqueComb(u,:);
        setvec = setvec(setvec>0);
        str = sprintf('%d→', setvec);
        str = str(1:end-1); % remove last arrow
        allStimSetStrings{end+1} = str;
    end
end

%% ===== Unique stimulation sets across all PTD folders =====
allStimSetStrings = unique(allStimSetStrings);

%% Depth map
d = Depth_s(Electrode_Type);
chRec = d(target_channel);

%%  ------------ PLOT ------------

for ss = 1:numel(allStimSetStrings)

    stimStr = allStimSetStrings{ss};
    stimVec = sscanf(strrep(stimStr,'→',' '), '%d')';

    figure('Color','w','Position',[200 50 950 1400]);
    tl = tiledlayout(numel(seq_folders)+1,1,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Sequential Stimulation — Ch %d — Set %s', ...
          target_channel, stimStr), 'FontSize',16);

    %% ====================== RASTERS ==========================
    for f = 1:numel(seq_folders)
        matchRow = [];
        for r = 1:Seq(f).nSets
            v = Seq(f).uniqueComb(r,:);
            v = v(v>0);
            if isequal(v,stimVec)
                matchRow = r;
                break;
            end
        end
        if isempty(matchRow)
            warning('Folder %d has no matching stimulation set for %s',f,stimStr);
            continue;
        end
        ax = nexttile(tl);
        hold(ax,'on'); box(ax,'off');
        sp_ch = Seq(f).sp{chRec};
        trig  = Seq(f).trig;
        PTDval= Seq(f).PTD(1);
        colorF= cmapPTD(f,:);
        trial_ids = find(Seq(f).combClass==matchRow & Seq(f).trialAmps==plot_amps);
        title(ax, sprintf('PTD = %g ms — Set %s', PTDval, stimStr));
        y = 0;
        for tr = trial_ids'
            t0 = trig(tr)/FS*1000;
            tt = sp_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;

            scatter(ax, tt, y*ones(size(tt)), 12, colorF, 'filled');
            y = y + 1;
        end
        xline(ax,0,'r--');
        xlim(ax,ras_win);
        ylim(ax,[0 max(y,1)]);
        ylabel(ax,'Trials');
    end

    %% ======================= PSTH ===========================
    ax = nexttile(tl);
    hold(ax,'on'); box(ax,'off');
    title(ax,'PSTH');

    for f = 1:numel(seq_folders)

        % find matching set row in this folder
        matchRow = [];
        for r = 1:Seq(f).nSets
            v = Seq(f).uniqueComb(r,:);
            v = v(v>0);
            if isequal(v,stimVec)
                matchRow = r;
                break;
            end
        end
        if isempty(matchRow), continue; end

        sp_ch = Seq(f).sp{chRec};
        trig  = Seq(f).trig;
        colorF= cmapPTD(f,:);
        PTDval= Seq(f).PTD(1);

        trial_ids = find(Seq(f).combClass==matchRow & Seq(f).trialAmps==plot_amps);
        if isempty(trial_ids), continue; end

        counts = zeros(1,length(edges)-1);

        for tr = trial_ids'
            t0 = trig(tr)/FS*1000;
            tt = sp_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            counts = counts + histcounts(tt,edges);
        end

        rate = filter(g,1, counts/(numel(trial_ids)*bin_s));

        plot(ax, ctrs, rate, 'LineWidth',2, 'Color', colorF, ...
             'DisplayName', sprintf('PTD %g ms', PTDval));
    end

    xline(ax,0,'r--','LineWidth',1,'HandleVisibility','off');
    xlabel(ax,'Time (ms)');
    ylabel(ax,'Rate (sp/s)');
    xlim(ax,ras_win);
    legend(ax,'Box','off');
end