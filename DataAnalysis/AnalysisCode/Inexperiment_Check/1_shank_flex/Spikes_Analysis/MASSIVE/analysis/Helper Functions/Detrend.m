% %% New detrend method
% function wave = Detrend(wave)
% FS = 30000;
% tapers = [0.025,400];
% n = tapers(1);
% w = tapers(2);
% p = n*w;
% k = floor(2*p-1);
% tapers = [n,p,k];
% tapers(1) = tapers(1).*FS;
% tapers = dpsschk(tapers);
% filt = mtfilt(tapers, FS, 0);
% filt = single(filt./sum(filt));
% tmp = single(medfilt1(double(wave),45)');
% wave = conv(tmp(1,46:end),filt);
% end

% function e_wave = Detrend(wave)
% d_wave = wave;
% BIN = 50*30; % 1 msec lines (possible want to choose this based on freq of interest?
% XSTART = find(wave == max(wave),1);
% MIN = find(wave(XSTART:end) == min(wave(XSTART:end)),1) + XSTART;
% LPEAK = find(wave(MIN:end) == max(wave(MIN:end)),1) + MIN;
% XSTOP = length(wave);
% BP = XSTART:BIN:XSTOP;
% d_wave(XSTART:XSTOP) = detrend(wave(XSTART:XSTOP),'linear',BP);
% MIN = find(d_wave(LPEAK:end) == min(d_wave(LPEAK:end)),1) + LPEAK;
% XSTOP = find(d_wave(MIN:end) == max(d_wave(MIN:end)),1) + MIN;
% X = LPEAK:1:XSTOP;
% CF = polyfit(X,d_wave(LPEAK:XSTOP),3);
% Y = polyval(CF,X);
% e_wave = d_wave;
% e_wave(LPEAK:XSTOP) = e_wave(LPEAK:XSTOP) - Y;
% figure;
% subplot(1,3,1);
% plot(wave);
% subplot(1,3,2);
% plot(d_wave);
% subplot(1,3,3);
% plot(e_wave);
% end

% function wave = Detrend(wave)
% %% Initial variables
% WINDOW = [-30,30]; % Sets up a 1 msec window on either side of the PoI
% START = -WINDOW(1)+1;
% STOP = length(wave)-WINDOW(2);
% 
% for P = START:STOP
%     % For every single point in the submitted wave
%     MEAN = median(wave(P+WINDOW(1):P+WINDOW(2)));
%     wave(P) = wave(P) - MEAN;
% end
% 
% end

%% Depreciated
% %% Function analyses a provided waveform and detrends post-artefact
function detrendedwave = Detrend(wave,ORDER)
% TRIG is the trigger line associated with this waveform
% NOTE - This may easily run into issues if given a waveform with more than
% one artefact in it. Be careful.
Y = zeros(1,length(wave));
XSTART = find(abs(wave(495*30:505*30)) == max(abs(wave(495*30:505*30))),1,'last') + 495*30;
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
detrendedwave = wave - Y;

% XSTART = find(temp_wave == max(temp_wave(505*30:end)),1);
% XSTOP = length(temp_wave);
% tmp = temp_wave(XSTART:XSTOP);
% % Fit a polynomial
% X = XSTART:1:XSTOP;
% [CF, S, mu] = polyfit(X,tmp,ORDER);
% Y(XSTART:XSTOP) = polyval(CF,X,S,mu);
% detrendedwave = temp_wave - Y;
end