%% Initial Setup
filepath = pwd;
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
if nChn>32
    E_Mapnumber=1;
else
    E_Mapnumber=0;
end
name = pwd;
[FP,name,ext] = fileparts(name);
name = name(1:end-14);

%% Deal with the digital lines
% Detect the falling edge of each trigger 
disp('Cleaning the digital lines');
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);
stimDig = flip(find(digital_in == 1)); % Fix for finding the falling line instead of the rising line
visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
dt = diff(visDig);
kill = dt == -1;
visDig(kill) = [];
nStamps = max([length(stimDig); length(visDig)]);
time_stamps = nan(1,nStamps);
time_stamps(1,1:length(stimDig)) = flip(stimDig);

%% Correct for 200 us delay
if ~isempty(dir('*exp_datafile*.mat'))
    Tname = dir('*exp_datafile*.mat');
    Tname = Tname.name;
    StimParams = load(Tname,'StimParams');
    StimParams = StimParams.StimParams;
    delay = StimParams(2,19);
    delay = delay{1};
else
    delay = 0;
end
time_stamps = time_stamps + delay*(1e-6)*FS;

%% Correct for jitter
d = Depth(E_Mapnumber); 
SKIP = 1;
try 
    loadStimChn;
catch 
    stimChn = 16;
end
trig = time_stamps;
if (SKIP == 0)
warning('off','signal:findpeaks:largeMinPeakHeight');
    for t = 1:length(trig)
        v = loadRAW(trig(t));
        [~,adj] = findpeaks(abs(diff(v(d(stimChn(1)),:))),'MinPeakHeight',500);
        if isempty(adj)
            continue;
        end
        % adj = adj(1) - (500*FS/1e3);
        adj = adj(1);
        trig(t) = trig(t) + adj;
        v = loadRAW(trig(t));
    end
end

trig_ms = (double(trig) ./ FS) * 1000;  % trigger time (ms)
trig_diff_ms = diff(trig_ms);
trig_diff_s = trig_diff_ms/1000;
mean_gap_ms = mean(trig_diff_ms);

% ---- Print trigger summary ----
fprintf('\n=== TRIGGER SUMMARY ===\n');
fprintf('Number of triggers: %d\n', length(trig));
fprintf('First trigger time: %.3f ms, %.3f s\n', trig_ms(1), trig_ms(1)/1000);
fprintf('Last  trigger time: %.3f ms, %.3f s, %.3f mins\n', trig_ms(end), trig_ms(end)/1000, trig_ms(end)/1000/60 );
if ~isempty(trig_diff_ms)
    fprintf('Inter-trigger gap (ms): mean=%.3f ms\n', mean_gap_ms);
else
    fprintf('Only one trigger present; no gaps to report.\n');
end

% plot trigger time difference 
figure();
plot(trig_diff_s);
xlabel('Trigger Number')
ylabel('Time (s)')


% trig_fid = fopen([name '.trig.dat'],'w');
% fwrite(trig_fid,trig,'double');
% fclose(trig_fid);


