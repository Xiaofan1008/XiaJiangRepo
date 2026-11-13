function [distance_um,note] = ChnPairDistance(chA,chB,Electrode_Type)
%CHNPAIRDISTANCE Physical distance (um) between two map-position channels.
%
% Uses ChnMap2Coord.m for each channel's (x,y) coordinate. Returns NaN,
% with an explanatory note, whenever the distance is not well-defined:
%   - either channel is an invalid map position for this Electrode_Type
%   - the two channels belong to different probe_group values (e.g. the
%     four-shank vs. single-shank probes under the 64chn+32chn hybrid,
%     Electrode_Type 3) - the physical offset between separate probes is
%     not known, so a Euclidean distance across them would be fabricated.
%
% OUTPUT
%   distance_um : scalar, NaN if undefined
%   note        : '' if distance_um is valid, otherwise a short reason

[xA,yA,groupA] = ChnMap2Coord(chA,Electrode_Type);
[xB,yB,groupB] = ChnMap2Coord(chB,Electrode_Type);

distance_um = NaN;
note = '';

if isempty(groupA) || isempty(groupB)
    note = 'invalid map position';
    return;
end

if ~strcmp(groupA,groupB)
    note = 'different probes - offset unknown';
    return;
end

if isnan(xA) || isnan(xB)
    note = 'x coordinate unknown for this probe';
    return;
end

distance_um = sqrt((xA-xB)^2+(yA-yB)^2);

end
