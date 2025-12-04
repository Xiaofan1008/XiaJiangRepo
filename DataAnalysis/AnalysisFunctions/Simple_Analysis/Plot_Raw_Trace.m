 %% Raw Trace Viewer (choose channels, amplitudes, PTDs, sets)
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% User settings
data_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim1_251202_121651';

channels_to_plot = [5 12 20];             % channels to plot (Intan numbering)
amps_to_plot     = [5 6 10];               % amplitudes to include (µA)
ptd_to_plot      = [];                    % PTDs to include (µs), empty = all
sets_to_plot     = [];                    % stimulation sets to include, empty = all

nTrials_to_plot  = 20;                    % number of trials to include
plot_window_ms   = [-5 40];               % window around trigger (ms)

Electrode_Type = 2;                       % 0 rigid, 1 flex, 2 4-shank flex

%% Check folder
if ~isfolder(data_folder), error('Invalid folder'); end
cd(data_folder);

%% Base name
parts = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');
if numel(u)>=4, base_name = lastfld(1:u(end-1)-1);
else, base_name = lastfld;
end

%% Intan info
fileinfo = dir('info.rhs');
[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);

%% Load amplifier.dat
fid = fopen('amplifier.dat','r');
if fid < 0, error('Cannot open amplifier.dat'); end

%% Load trigger file
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

%% Load StimParams
fDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fDIR),'No *_exp_datafile_*.mat found.');
S = load(fDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
n_Trials = S.n_Trials;
E_MAP = S.E_MAP;

%% Amplitudes
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps,~,ampIdx] = unique(trialAmps(:));
Amps(Amps==-1) = 0;
n_AMP = numel(Amps);

%% PTD decode
postTrig_all = cell2mat(StimParams(2:end,6));
postTrig = postTrig_all(2:simultaneous_stim:end);
[PTDs,~,ptdIdx] = unique(postTrig);
n_PTD = numel(PTDs);

%% Sequence-sensitive stimulation sets
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = idx_all(rr);
    v = v(v>0);
    stimChPerTrial_all{t} = v(:)';
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t, 1:numel(v)) = v;
end

[uniqueComb,~,combClass] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);
combClass_win = combClass;

%% PTD amplitude set selection
if isempty(ptd_to_plot), ptd_sel = PTDs;
else, ptd_sel = intersect(PTDs, ptd_to_plot);
end
if isempty(sets_to_plot), set_sel = 1:nSets;
else, set_sel = sets_to_plot;
end
if isempty(amps_to_plot), amps_sel = Amps;
else, amps_sel = intersect(Amps, amps_to_plot);
end

%% Time window and buffer
samp_win = round(plot_window_ms/1000 * FS);
Nsamp = samp_win(2)-samp_win(1)+1;
time_ms = (samp_win(1):samp_win(2)) / FS * 1000;

%% Channel mapping
d = Depth_s(Electrode_Type);

%% Main loop
for ch_plot = channels_to_plot

    ch_intan = d(ch_plot);

    for set_id = set_sel

        stimVec = uniqueComb(set_id,:);
        stimVec = stimVec(stimVec>0);
        set_label = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimVec,'UniformOutput',false),'→');

        for ptd_val = ptd_sel
            tidx_ptd = find(postTrig == ptd_val);

            for amp_val = amps_sel
                tidx_amp = find(trialAmps == amp_val);

                tlist = intersect(tidx_ptd, tidx_amp);
                tlist = intersect(tlist, find(combClass_win == set_id));

                if isempty(tlist), continue; end

                tlist = tlist(1:min(nTrials_to_plot,numel(tlist)));

                figName = sprintf('RawTrace_Ch%d | Set %s | Amp %d | PTD %dus', ...
                                    ch_plot, set_label, amp_val, ptd_val);
                figure('Name',figName,'Color','w','Position',[200 200 1200 600]);
                hold on

                for k = 1:numel(tlist)
                    tr = tlist(k);
                    start_idx = trig(tr) + samp_win(1);
                    byte_pos = start_idx * nChn * 2;

                    fseek(fid, byte_pos, 'bof');
                    raw = fread(fid, [nChn, Nsamp], 'int16') * 0.195;

                    plot(time_ms, raw(ch_intan,:), 'LineWidth', 1.0);
                end

                xline(0,'r--');
                xlabel('Time (ms)');
                ylabel('µV');
                title(sprintf('Raw trace | Ch %d | %s | %d µA | PTD %d µs', ...
                       ch_plot, set_label, amp_val, ptd_val));
                grid on
            end
        end
    end
end

fclose(fid);