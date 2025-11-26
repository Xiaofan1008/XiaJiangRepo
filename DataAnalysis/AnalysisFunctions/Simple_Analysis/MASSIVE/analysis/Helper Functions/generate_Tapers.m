%% General function for generating filters for use in extracting LFP and MUA activity from raw neural data
function [filt] = generate_Tapers(FS,len,BW,CF)
%% Instructions
% Takes two parameters and returns a filter. _len represents the sequence
% length in seconds, and _BW represents the effective bandwidth.
% FS is the sampling frequency of the raw datafile.
% Requires the presence of dpsschk and mtfilt in MATLAB's path.

% tapers = [sequence length (seconds), effective bandwidth (Hz)]
%% Initialisation of parameters
if (nargin < 4)
    FS = 30000;
    len = 0.0025;
    BW = 400;
    CF = 0;
end

%% FILTER
tapers = [len,BW];
n = tapers(1);
w = tapers(2);
p = n*w;
k = floor(2*p-1);
tapers = [n,p,k];
tapers(1) = tapers(1).*FS;
tapers = dpsschk(tapers);
filt = mtfilt(tapers, FS, CF);
filt = single(filt./sum(filt)); % Converts output to be non-complex

end