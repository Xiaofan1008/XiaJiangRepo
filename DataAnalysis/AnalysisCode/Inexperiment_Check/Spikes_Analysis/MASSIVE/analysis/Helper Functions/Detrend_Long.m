%% Function analyses a provided waveform and detrends post-artefact
function detrendedwave = Detrend_Long(wave,ORDER)
% TRIG is the trigger line associated with this waveform
% NOTE - This may easily run into issues if given a waveform with more than
% one artefact in it. Be careful.
Y = ones(1,length(wave))*mean(wave(10:90));
%Y = zeros(1,length(wave));
XSTART = find(abs(wave(149*30:160*30)) == max(abs(wave(149*30:154*30))),1,'last') + 148*30;
XSTOP = length(wave);
% Downsample wave
try
    tmp = wave(XSTART:XSTOP);
catch
    disp('Help!');
end
% Fit a polynomial
X = XSTART:1:XSTOP;
[CF, S, mu] = polyfit(X,tmp,ORDER);
Y(XSTART:XSTOP) = polyval(CF,X,S,mu);
detrendedwave = wave;
detrendedwave(XSTART:XSTOP) = wave(XSTART:XSTOP) - Y(XSTART:XSTOP);
end