%% Raw / Denoised / Filtered Trace Viewer (sequence-sensitive sets, separate figs for set×amp, subplots for PTDs)
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

channels_to_plot = 1:32;                % channels to plot (Depth_s index)
amps_to_plot     = [5,6];                 % amplitudes to include (µA)
ptd_to_plot      = [];                  % PTDs (ms), [] means all
sets_to_plot     = [];                  % stimulation sets, [] means all

nTrials_to_plot  = 5;                  % how many trials to draw per condition
plot_window_ms   = [-5 80];             % window around trigger

Electrode_Type   = 1;                   % 0 rigid, 1 flex, 2 4-shank flex

% 'raw'  = amplifier.dat
% 'dn'   = amplifier_dn_sab.dat   (denoised)
% 'mu'   = <base_name>.mu_sab.dat (filtered / MUA)
trace_type = 'raw';    % raw trace
% trace_type = 'dn';    % artifact blanked trace
% trace_type = 'mu';    % bandpass filtered trace

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

        for i_amp = 1:length(amps_sel)
            amp_val = amps_sel(i_amp);

            % ---- figure name + title depend on trace type ----
            figName = sprintf('%sTrace_Ch%d_Set%s_Amp%d', ...
                              data_label, ch_plot, set_label, amp_val);
            figure('Name',figName,'Color','w','Position',[150 150 1400 800]);

            nPTD  = length(ptd_sel);
            nRows = ceil(sqrt(nPTD));
            nCols = ceil(nPTD / nRows);

            for i_ptd = 1:nPTD
                ptd_val = ptd_sel(i_ptd);

                tidx_set = find(combClass_win == set_id);
                tidx_ptd = find(ismember(postTrig, ptd_val));
                tidx_amp = find(trialAmps == amp_val);

                tlist = intersect(intersect(tidx_set, tidx_ptd), tidx_amp);

                if isempty(tlist), continue; end
                tlist = tlist(1:min(nTrials_to_plot, numel(tlist)));

                subplot(nRows, nCols, i_ptd); hold on

                for k = 1:numel(tlist)
                    tr = tlist(k);
                    start_idx = trig(tr) + samp_win(1);
                    byte_pos  = start_idx * nChn * 2;

                    fseek(fid, byte_pos, 'bof');
                    data_block = fread(fid, [nChn, Nsamp], 'int16') * 0.195;

                    plot(time_ms, data_block(ch_intan,:), 'LineWidth', 1);
                end

                xline(0,'r--');
                title(sprintf('Inter-Stimulus-Interval %d ms (%d trials)', ...
                      ptd_val/1000, numel(tlist)));
                xlabel('Time (ms)');
                ylabel('µV');
                ylim([-7000,7000]);
            end

            sgtitle(sprintf('%s Trace | Ch %d | Set %s | %d µA', ...
                data_label, ch_plot, set_label, amp_val));
        end
    end
end

fclose(fid);