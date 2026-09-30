%% Check 2: amplifier.dat vs Strobe_new.mu_sab.dat, plus 50 Hz line-noise check
clear; close all;

% ---------------- settings ----------------
folder = '/Volumes/MACData/Data/Data_Xia/DX012/Strobe_new_251125_211617';
fs   = 30000;
nCh  = 64;
gain = 0.195;    % µV per bit
row  = 10;       % channel row in the file to inspect (try a few)
t0   = 60;       % start time of segment (s)
dur  = 20;       % segment length (s)

% ---------------- file sizes ----------------
a  = dir(fullfile(folder, 'amplifier.dat'));
m  = dir(fullfile(folder, 'Strobe_new.mu_sab.dat'));
tf = dir(fullfile(folder, 'time.dat'));

nAmp = a.bytes / (2*nCh);
fprintf('amplifier.dat : %.0f samples = %.1f min\n', nAmp, nAmp/fs/60);
fprintf('time.dat      : %.0f samples (should match)\n', tf.bytes/4);
bpf = m.bytes / nAmp;
fprintf('mu_sab.dat    : %.3f bytes per amplifier sample\n', bpf);
fprintf('                (128 = 64 ch int16, 256 = 64 ch single; other = different format)\n');

if nAmp/fs < t0 + dur
    t0 = 0;  dur = min(dur, floor(nAmp/fs));
end

% ---------------- read segments ----------------
amp = read_seg(fullfile(folder, 'amplifier.dat'), nCh, fs, t0, dur, 'int16') * gain;

if abs(bpf - 2*nCh) < 1e-6
    mu = read_seg(fullfile(folder, 'Strobe_new.mu_sab.dat'), nCh, fs, t0, dur, 'int16');
elseif abs(bpf - 4*nCh) < 1e-6
    mu = read_seg(fullfile(folder, 'Strobe_new.mu_sab.dat'), nCh, fs, t0, dur, 'single');
else
    mu = [];
    warning('mu_sab.dat layout not recognised - skipping it. Send Claude the bytes-per-sample number.');
end

xa = amp(row,:) - mean(amp(row,:));
t  = (0:numel(xa)-1) / fs;
one = 1:fs;   % first second for trace plots

% ---------------- spectra ----------------
win = 2*fs;                                   % 2 s windows -> ~0.5 Hz resolution
[pa, f] = pwelch(xa, win, win/2, [], fs);
if ~isempty(mu)
    xm = mu(row,:) - mean(mu(row,:));
    [pm, ~] = pwelch(xm, win, win/2, [], fs);
end

% ---------------- 50 Hz line-noise report ----------------
fprintf('\nLine-noise peaks in amplifier.dat, row %d (peak / neighbouring power):\n', row);
for f0 = 50:50:300
    pk = max(pa(abs(f - f0) < 0.6));
    nb = median(pa(abs(f - f0) > 2 & abs(f - f0) < 6));
    fprintf('  %3d Hz : %6.1f x\n', f0, pk/nb);
end
fprintf('  (~1-2 x = no line noise; > ~5 x = clear line noise)\n');

% ---------------- plots ----------------
figure('Color', 'w', 'Position', [100 80 1100 850]);

subplot(4,1,1);
plot(t(one), xa(one), 'k');
title(sprintf('amplifier.dat, row %d (1 s)', row)); ylabel('µV'); xlabel('s');

subplot(4,1,2);
if ~isempty(mu)
    plot(t(one), xm(one), 'b');
    title(sprintf('mu\\_sab.dat, row %d (1 s)', row)); ylabel('raw units'); xlabel('s');
else
    text(0.5, 0.5, 'mu\_sab.dat not read (unknown format)', 'HorizontalAlignment', 'center');
    axis off;
end

subplot(4,1,3);
loglog(f, pa/max(pa), 'k'); hold on;
if ~isempty(mu), loglog(f, pm/max(pm), 'b'); legend('amplifier', 'mu\_sab'); end
xline(300, 'r:'); xlim([0.5 10000]); grid on;
title('Power spectrum (normalised)'); xlabel('Hz'); ylabel('power');

subplot(4,1,4);
sel = f >= 20 & f <= 320;
semilogy(f(sel), pa(sel), 'k'); hold on;
for f0 = 50:50:300, xline(f0, 'r:'); end
grid on; xlim([20 320]);
title('amplifier.dat zoom: red lines = 50 Hz and harmonics'); xlabel('Hz'); ylabel('µV^2/Hz');

% ---------------- local function (must stay at end of file) ----------------
function x = read_seg(file, nCh, fs, t0, dur, prec)
    switch prec
        case 'int16',  bps = 2;
        case 'single', bps = 4;
    end
    fid = fopen(file, 'r');
    assert(fid > 0, 'Cannot open %s', file);
    fseek(fid, round(t0*fs) * nCh * bps, 'bof');
    x = fread(fid, [nCh, round(dur*fs)], [prec '=>double']);
    fclose(fid);
end