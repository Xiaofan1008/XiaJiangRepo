function [x_um,y_um,probe_group] = ChnMap2Coord(map_position,Electrode_Type)
%CHNMAP2COORD Physical (x,y) coordinate of a map-position channel, in um.
%
% GEOMETRY (fixed by probe design, confirmed against the probe
% schematics):
%   Single-shank probe (32 ch): one column, 50 um pitch, channel 1 = tip
%     (deepest), channel 32 = most superficial.
%   Four-shank probe (64 ch): four columns of 16 ch each, 50 um pitch
%     within a shank, 200 um between adjacent shanks (x = 0/200/400/600
%     um for shank 1/2/3/4, in that physical left-to-right order).
%     Channel 1 of each shank's range = tip (deepest).
%       Shank 1: map positions 1-16
%       Shank 2: map positions 33-48
%       Shank 3: map positions 49-64
%       Shank 4: map positions 17-32
%
% Electrode_Type selects which geometry profile applies (same convention
% as ChnMap.m):
%   0, 1 : single-shank only (32 ch)
%   2    : four-shank only (64 ch)
%   3    : 64chn+32chn hybrid - map positions 1-64 are the four-shank
%          probe, 65-96 are the single-shank probe. These are TWO
%          SEPARATE probes (different ports/insertion sites); the
%          physical offset between them is not known, so x_um/y_um for
%          the single-shank group under this type are returned NaN, and
%          probe_group lets callers refuse to compute a cross-probe
%          distance rather than silently returning a wrong number.
%
% OUTPUT
%   x_um, y_um   : physical coordinates in micrometres, NaN if
%                  map_position is invalid for this Electrode_Type, or if
%                  the coordinate is not known (hybrid single-shank x).
%   probe_group  : 'four_shank' or 'single_shank' (string), '' if
%                  map_position is invalid.

x_um = NaN;
y_um = NaN;
probe_group = '';

map_position = double(map_position);
if ~isscalar(map_position) || ~isfinite(map_position) || map_position < 1 || ...
        fix(map_position) ~= map_position
    return;
end

SHANK_PITCH_UM = 200;
ROW_PITCH_UM   = 50;

switch Electrode_Type
    case {0,1}   % single-shank only, 32 ch
        if map_position > 32
            return;
        end
        probe_group = 'single_shank';
        x_um = 0;
        y_um = (map_position-1)*ROW_PITCH_UM;

    case 2       % four-shank only, 64 ch
        [x_um,y_um] = four_shank_coord(map_position,SHANK_PITCH_UM,ROW_PITCH_UM);
        if ~isnan(x_um)
            probe_group = 'four_shank';
        end

    case 3       % 64chn+32chn hybrid
        if map_position <= 64
            [x_um,y_um] = four_shank_coord(map_position,SHANK_PITCH_UM,ROW_PITCH_UM);
            if ~isnan(x_um)
                probe_group = 'four_shank';
            end
        elseif map_position <= 96
            probe_group = 'single_shank';
            row = map_position-64;
            x_um = NaN;   % offset between the two probes is not known
            y_um = (row-1)*ROW_PITCH_UM;
        end

    otherwise
        error('ChnMap2Coord: unsupported Electrode_Type %d.',Electrode_Type);
end

end

function [x_um,y_um] = four_shank_coord(map_position,shank_pitch_um,row_pitch_um)

x_um = NaN;
y_um = NaN;

if map_position >= 1 && map_position <= 16
    shank_index = 1;
    row = map_position;
elseif map_position >= 33 && map_position <= 48
    shank_index = 2;
    row = map_position-32;
elseif map_position >= 49 && map_position <= 64
    shank_index = 3;
    row = map_position-48;
elseif map_position >= 17 && map_position <= 32
    shank_index = 4;
    row = map_position-16;
else
    return;
end

x_um = (shank_index-1)*shank_pitch_um;
y_um = (row-1)*row_pitch_um;

end
