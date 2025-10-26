function depth = Depth_s(caseType)
% this function map the electrode ID to the file channels
%   caseType - 0: single shank rigid, p=3
%              1: single shank flex,  p=5
%              2: four shank flex,    p=6

% Read amplifier channel info
amplifier_channels = read_Intan_RHS2000_file;
nChn = length(amplifier_channels);

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
    otherwise
        error('Invalid caseType. Use 0 (rigid), 1 (flex), or 2 (4-shank flex).');
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



