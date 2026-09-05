function Snips_clean_uV = blank_Artifacts_fixed(pre_ms, post_ms, inFolder, inFile) 
%   Interpolate (blank) stimulation artifacts in amplifier.dat around each trigger
%
%   Usage:
%   blankArtifactsToFile();                          % defaults: [-1, +2] ms, pwd, amplifier.dat -> amplifier_blank.dat
%   blankArtifactsToFile(1, 2);                      % set window
%   blankArtifactsToFile(1, 2, 'D:\data\rec1');      % custom folder
%   blankArtifactsToFile(1, 2, [], 'amplifier.dat'); % custom names

% ------------ defaults ------------
if nargin < 1 || isempty(pre_ms),  pre_ms  = 1;   end    % ms before trigger
if nargin < 2 || isempty(post_ms), post_ms = 2;   end    % ms after trigger
if nargin < 3 || isempty(inFolder), inFolder = pwd; end
if nargin < 4 || isempty(inFile),   inFile   = 'amplifier.dat'; end

inPath  = fullfile(inFolder, inFile);

assert(exist(inPath,'file')==2, 'Input file not found: %s', inPath);

% ---- Intan information ---- 
[amplifier_channels, frequency_parameters] = read_Intan_RHS2000_file;
nChn = numel(amplifier_channels);
FS   = frequency_parameters.amplifier_sample_rate;

% ----- triggers ----
trig = loadTrig(0); %  (sample indices)
nTrig = length(trig);
assert(~isempty(trig), 'No triggers found.');

% ---- file dimensions ----
info = dir(inPath);
BytePerSample = 2;                           % bytes per int16 sample
total_samples_per_ch = info.bytes / (nChn * BytePerSample);
total_samples_per_ch = floor(total_samples_per_ch);

% ---- data read window: -10 ms +100 ms ----
data_pre_ms  = 10;  
data_post_ms = 100;
data_preSamp = round(data_pre_ms  / 1000 * FS);
data_postSamp = round(data_post_ms / 1000 * FS);
data_len = data_preSamp + data_postSamp + 1;          
t_ms = (-data_preSamp:data_postSamp) / FS * 1000; 

 % ---- Artifact window ----
art_preSamp  = round(pre_ms / 1000 * FS);
art_postSamp = round(post_ms / 1000 * FS);

% --- Preallocate ---
Snips_clean_uV = [];
usedTrig = [];
tmpCell = cell(numel(trig), 1);
keepCount = 0;
fid = fopen(inFile, 'r');
assert(fid ~= -1, 'Cannot open file: %s', inFile);

for t = 1:nTrig
    thisTrig = trig(t);
    left  = thisTrig - data_preSamp;
    right = thisTrig + data_postSamp;

    % Skip triggers near edges
    if left < 1 || right > total_samples_per_ch
        continue;
    end

    nRead = data_len;
    fseek(fid, int64(left-1) * nChn * bytes_per_sample, 'bof');
    seg = fread(fid, [nChn, nRead], 'int16') * 0.195;  % µV

    % Artifact indices inside seg
    art_start = data_preS - art_preSamp + 1;  % relative to seg start
    art_end   = data_preS + art_postSamp + 1;

    for ch = 1:nChn
        yL = seg(ch, art_start-1);  
        yR = seg(ch, art_end+1);  
        ramp = linspace(0, 1, art_end - art_start + 1);
        seg(ch, art_start:art_end) = yL + (yR - yL) .* ramp;
    end

    keepCount = keepCount + 1;
    tmpCell{keepCount} = seg; 
    usedTrig(keepCount,1) = ti;
end

if keepCount > 0
    Snips_clean_uV = cat(3, tmpCell{1:keepCount}); % [nChn x big_len x nTrig]
else
    Snips_clean_uV = zeros(nChn, big_len, 0, 'double');
end
end