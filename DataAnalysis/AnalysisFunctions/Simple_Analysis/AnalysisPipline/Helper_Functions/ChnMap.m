function depth = ChnMap(caseType, nChn)
% ChnMap  Map electrode ID to file channels, WITHOUT needing the raw
%         Intan recording files on disk.
%
%   depth = ChnMap(caseType, nChn)
%
%   Identical to Depth_s(caseType), except the channel count is passed
%   in directly (nChn) instead of being read from amplifier.dat via
%   read_Intan_RHS2000_file. The actual electrode-to-channel mapping
%   itself only ever came from ProbeMAP/E_MAP (fixed, deterministic --
%   has nothing to do with any specific recording's data), so this is
%   the only thing that needed to change. Use this instead of Depth_s
%   whenever the raw Intan files for a dataset may no longer exist
%   (e.g. amplifier.dat deleted to save space) but you already know the
%   channel count from elsewhere (e.g. numel(sp_waveforms) or
%   numel(sp_clipped) from that dataset's own spike files).
%
%   caseType - 0: single shank rigid, p=3
%              1: single shank flex,  p=5
%              2: four shank flex,    p=6
%              3: 64chn+32chn hybrid, p=7
%   nChn     - number of amplifier channels (NOT read from any file --
%              you supply this directly, e.g. from numel(sp_waveforms)).

% Load mapping
E_MAP = ProbeMAP;

% Select map index p based on caseType
switch caseType
    case 0
        p = 3; % single shank rigid
    case 1
        p = 5; % single shank flex
    case 2
        p = 6; % four shank flex
    case 3
        p = 7; % 64chn+32chn hybrid
    otherwise
        error('Invalid caseType. Use 0 (rigid), 1 (flex), 2 (4-shank flex), or 3 (hybrid).');
end

% Map each channel to a depth
depth = zeros(1, nChn);
for n = 2:nChn+1
    str = E_MAP{n, p};
    if str(1) == 'A'
        str = str2double(str(4:5));
    elseif str(1) == 'B'
        str = str2double(str(4:5)) + 32;
    elseif str(1) == 'C'
        str = str2double(str(4:5)) + 64;
    elseif str(1) == 'D'
        str = str2double(str(4:5)) + 96;
    else
        str = 0;
    end
    depth(n-1) = str;
end
if any(depth == 0)
    depth = depth + 1;
end
depth = depth';  % return as column vector
end